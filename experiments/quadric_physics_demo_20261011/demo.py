"""Archived quadric Heart: material visualization and measured Bullet replay.

The voxel panel samples the SAME final quadric; it is never a voxel optimization.
Only the original outer hull collides. Rejected floating interior derivatives
are display-only; physical mass properties come from rational saved material.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import subprocess
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'experiments/faster_heart_contact_20261011'),
             str(ROOT/'experiments/heart_stem_contact')]
from contact_matrix import physical_properties, scenarios

SOURCE=ROOT/'.superpowers/sdd/2026-10-07-spin2014-objectives-quadric/heart-lowering-seed1-normalized-replay-run-01/output'

def write(path,data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False))

def authenticated_bytes(folder,record):
    path=folder/record['name']
    if path.parent!=folder or path.is_symlink():raise ValueError('INVALID_PAYLOAD_PATH')
    data=path.read_bytes()
    if len(data)!=record['size'] or hashlib.sha256(data).hexdigest()!=record['sha256']:
        raise ValueError('SAVED_PAYLOAD_CHANGED')
    return data

def to_physics(moments):
    return physical_properties(moments,mass_unit_kg=.01875,length_unit_m=.025)

def prepare(out):
    import io
    import trimesh
    from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
    folder=SOURCE/'material-geometry';manifest=json.loads((folder/'manifest.json').read_text())
    audit=json.loads((SOURCE/'independent-material-audit.json').read_text())['audit']
    run=json.loads((SOURCE/'native-run.json').read_text())
    comparison=json.loads((SOURCE.parent.parent/'heart-lowering-seed1-normalized-replay-comparison-01/comparison.json').read_text())
    assert comparison['identical_material_payloads'] and comparison['identical_scientific_events']
    assert run['solve']['status']=='CONVERGED' and audit['independent_raw_residual']<1e-12
    assert audit['fixed_density']==audit['design_density']==1
    meshes={};raw={};hashes={}
    for name,description in manifest['meshes'].items():
        payload={}
        for field in ('vertices','faces','rational_vertices'):
            r=description.get(field)
            if r:
                data=authenticated_bytes(folder,r);hashes[str(folder/r['name'])]=hashlib.sha256(data).hexdigest()
                payload[field]=json.loads(data) if field=='rational_vertices' else np.load(io.BytesIO(data),allow_pickle=False)
        meshes[name]=(payload['vertices'],payload['faces'])
        if name in ('outer','material'):
            print('Independent exact integration:',name,flush=True)
            coords=payload.get('rational_vertices',payload['vertices'].tolist())
            values=integrate_exact_rational_raw_moments(coords,payload['faces'].tolist(),arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            raw[name]=[float(x) for x in values]
    np.testing.assert_allclose(raw['material'],audit['total_moments_model'],rtol=1e-12,atol=1e-14)
    v,f=meshes['outer'];outer=trimesh.Trimesh(v*.025,f,process=False)
    if not outer.is_watertight or not outer.is_winding_consistent or outer.volume<=0:raise ValueError('INVALID_OUTER')
    outer.export(out/'outer-si.obj',digits=17)
    hull=outer.convex_hull
    if not hull.is_watertight or hull.volume<=0:raise ValueError('INVALID_COLLISION_HULL')
    hull.export(out/'collision-hull-si.obj',digits=17)
    # Convex hull of source vertices preserves every plane support value.
    from scipy.spatial.transform import Rotation
    support=[]
    for angle in (0,2,5,15,45,90,135,180):
        for az in range(0,360,45):
            R=Rotation.from_rotvec(np.deg2rad(angle)*np.array([np.cos(np.deg2rad(az)),np.sin(np.deg2rad(az)),0])).as_matrix()
            support.append(abs((hull.vertices@R.T)[:,2].min()-(outer.vertices@R.T)[:,2].min()))
    if max(support)>1e-14:raise ValueError('OUTER_HULL_SUPPORT_MISMATCH')
    props={name:to_physics(raw[key]) for name,key in [('original','outer'),('optimized','material')]}
    info={name:dict(mass_kg=float(p[0]),center_m=p[1].tolist(),inertia_kg_m2=p[2].tolist(),
                   principal_moments=p[3].tolist(),principal_axes=p[4].tolist()) for name,p in props.items()}
    write(out/'material-properties.json',dict(materials=info,raw_moments=raw,source=str(SOURCE),
        source_sha256=hashes,independent_raw_residual=audit['independent_raw_residual'],
        accepted_updates=audit['accepted_updates'],objective='SPIN2014_TOP_LOWERING_ONLY',
        solver_method='FASTER_QUADRIC_SPHERICAL_NEWTON_CGAL',
        max_plane_support_error_m=max(support),collision='COMMON_OUTER_CONVEX_HULL',
        collision_scope='Plane support validated; contact manifold and real friction are not calibrated.',
        interior_float_mesh_use='DISPLAY_ONLY_NEVER_COLLISION',manufacturing_approved=False))
    return meshes,props,run,hashes

def figures(out,meshes,props,run):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.patches import Polygon
    from matplotlib.path import Path as PlotPath
    from matplotlib.collections import LineCollection
    import importlib.util
    spec=importlib.util.spec_from_file_location('saved_renderer',ROOT/'experiments/faster_quadric_visualization_20261011/render.py')
    renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)
    plt.rcParams['font.family']=font_manager.FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
    plt.rcParams['axes.unicode_minus']=False
    fig=plt.figure(figsize=(15,8),facecolor='#f4f7fc')
    gs=fig.add_gridspec(2,3,left=.045,right=.97,top=.88,bottom=.10,hspace=.4,wspace=.27)
    fig.suptitle('同一个 Heart：二次曲面如何改变内部质量分布',fontsize=22,x=.045,ha='left',color='#142944')
    for column,label in enumerate(('原始实心体','实际优化材料 · 展示性剖切')):
        ax=fig.add_subplot(gs[0,column],projection='3d')
        if column==0:
            v,f=meshes['outer'];tris=v[f];color='#dc8450'
            ax.add_collection3d(Poly3DCollection(tris,facecolor=color,edgecolor='none',alpha=.88))
        else:
            for name,color in [('fixed','#ccd5e1'),('retained','#16a6a1')]:
                ax.add_collection3d(Poly3DCollection(renderer.clip_triangles(*meshes[name],axis=1,value=.015),facecolor=color,edgecolor='none',alpha=.95))
        ax.set(xlim=(-1,1),ylim=(-.5,.5),zlim=(-.1,1.9),title=label)
        ax.set_box_aspect((2,1,2));ax.view_init(elev=15,azim=-70);ax.set_axis_off()
        p=props['original' if column==0 else 'optimized'];c=p[1]/.025
        ax.scatter(*c,color='#e23163',s=36);ax.plot([0,0],[0,0],[0,1.9],'--',c='#283e61')
    q=np.array(run['solve']['quadric_coefficients'])
    x=np.linspace(-.95,.95,401);z=np.linspace(.01,1.9,401);X,Z=np.meshgrid(x,z)
    F=renderer.polynomial(q,X,.017,Z)
    # Polygon bounds only determine display sampling, never simulation inertia.
    import trimesh
    v,f=meshes['design'];section=trimesh.Trimesh(v,f,process=False).section([0,1,0],[0,.017,0])
    loops=section.discrete
    def mask(a,b):
        points=np.column_stack([a.ravel(),b.ravel()]);inside=np.zeros(len(points),bool)
        for loop in loops:inside^=PlotPath(loop[:,[0,2]]).contains_points(points)
        return inside.reshape(a.shape)
    inside=mask(X,Z)
    ax=fig.add_subplot(gs[0,2]);ax.contourf(X,Z,np.where(inside,F,np.nan),levels=[-10,0,10],colors=['#f5d4aa','#a5e1d9'])
    ax.contour(X,Z,F,levels=[0],colors=['#168c88'],linewidths=2)
    for loop in loops:ax.plot(loop[:,0],loop[:,2],color='#5d718c')
    ax.set(title='Faster 表示：全局二次曲面',xlim=(-.95,.95),ylim=(0,1.9));ax.set_aspect('equal')
    ax.text(.02,.02,'10 个归一化系数 / 9 个独立方向',transform=ax.transAxes,fontsize=9)
    ax=fig.add_subplot(gs[1,0]);grid=np.linspace(-.95,.95,25);gz=np.linspace(0,1.9,25)
    cx,cz=np.meshgrid((grid[:-1]+grid[1:])/2,(gz[:-1]+gz[1:])/2)
    raster=np.where(mask(cx,cz),renderer.polynomial(q,cx,.017,cz)>0,np.nan)
    ax.pcolormesh(grid,gz,raster,cmap=matplotlib.colors.ListedColormap(['#f5d4aa','#a5e1d9']),edgecolors='#8293a8',linewidth=.4)
    ax.contour(X,Z,F,levels=[0],colors=['#137f7e'],linewidths=1.4)
    ax.set(title='体素表示示意：每个单元一个变量',xlim=(-.95,.95),ylim=(0,1.9));ax.set_aspect('equal')
    ax.text(.02,.02,'同一曲面的采样，未执行体素优化',transform=ax.transAxes,fontsize=9)
    trace=run['solve']['trace'];ax=fig.add_subplot(gs[1,1])
    ax.semilogy(range(len(trace)),[max(t['raw_equation_residual'],1e-18) for t in trace],color='#098f8b')
    ax.axhline(1e-12,ls='--',color='#da8150');ax.set(title='保存的真实 Newton 轨迹',xlabel='记录序号（含重启）',ylabel='原始方程残差');ax.grid(alpha=.2)
    ax=fig.add_subplot(gs[1,2]);ax.axis('off')
    before,after=props['original'],props['optimized']
    ax.text(0,.95,f'实心质量   {before[0]*1000:.3f} g\n优化质量   {after[0]*1000:.3f} g\n减重          {(1-after[0]/before[0])*100:.2f}%\n\n质心高度   {before[1][2]*1000:.3f} → {after[1][2]*1000:.3f} mm\n\n外形相同 · 内部材料不同\n粉点：实际质心\n灰色：固定壳，青色：保留材料',va='top',fontsize=12,color='#142944')
    fig.text(.045,.025,'重建方法案例 | Faster 二次曲面求解 + Spin-It 2014 惯量与平方高度目标 | 内部网格仅用于显示',fontsize=10)
    fig.savefig(out/'method-comparison.png',dpi=150);plt.close(fig)

def simulate(out,props):
    import trimesh
    from rounded_contact import simulate_scene
    protocol=json.loads((ROOT/'experiments/faster_heart_contact_20261011/rounded-asymmetric-v2/contact-protocol.json').read_text())
    # Existing backend validates the shared protocol; this delivery explicitly
    # runs ONE hull proxy and reports 32/32 of its own selected scenarios.
    selected=[x for x in scenarios(protocol) if x['proxy']=='default']
    write(out/'simulation-protocol.json',dict(physical_parameters=protocol,scenarios=selected,
        selection='32 scenes, common outer convex hull; no robust/VHACD comparison claimed',
        video_selection='azimuth 0 degrees, dt .00025, both materials; fixed before simulation'))
    outer=out/'outer-si.obj';vertices=trimesh.load(outer,force='mesh',process=False).vertices
    results=[]
    for i,scene in enumerate(selected):
        print(f'PyBullet {i+1}/{len(selected)} {scene["case_id"]}',flush=True)
        try:
            result=simulate_scene(outer,[out/'collision-hull-si.obj'],props[scene['material']],vertices,scene,protocol)
            result['proxy_description']='COMMON_ORIGINAL_OUTER_VERTEX_CONVEX_HULL'
            result['source_case']='ARCHIVED_HEART_SPIN2014_OBJECTIVE_QUADRIC'
        except Exception as error:
            result=dict(scene=scene,status='FAILED',error_type=type(error).__name__,error=str(error))
        write(out/(scene['case_id']+'.json'),result);results.append(result)
    summary=[]
    for r in results:
        summary.append({k:v for k,v in r.items() if k!='records'})
    write(out/'simulation-summary.json',dict(expected_scenarios=32,completed=sum(r['status']!='FAILED' for r in results),
        results=summary,scope='CONTROLLED_UNCALIBRATED_CONTACT_NOT_REAL_STOPPING_TIME_OR_MANUFACTURING'))
    return results

def render_video(out,props,results):
    import pybullet as p
    from bullet_contact_common_frame import create_body
    from PIL import Image,ImageDraw,ImageFont
    chosen=[next(r for r in results if r['scene']['material']==m and r['scene']['azimuth_degrees']==0 and r['scene']['dt_seconds']==.00025) for m in ('original','optimized')]
    if any(r['status']=='FAILED' for r in chosen):raise ValueError('PREREGISTERED_VIDEO_SCENE_FAILED')
    client=p.connect(p.DIRECT);bodies=[]
    p.setGravity(0,0,-9.81,physicsClientId=client)
    floor=p.createVisualShape(p.GEOM_BOX,halfExtents=[.4,.4,.001],rgbaColor=[.87,.9,.95,1],physicsClientId=client)
    p.createMultiBody(baseMass=0,baseVisualShapeIndex=floor,basePosition=[0,0,-.001],physicsClientId=client)
    colors=[(.91,.48,.25,1),(.02,.62,.60,1)]
    for m,color in zip(('original','optimized'),colors):
        bodies.append(create_body(p,client,out/'outer-si.obj',[out/'collision-hull-si.obj'],props[m],np.eye(3),np.zeros(3),35,color))
    fontpath='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
    font=ImageFont.truetype(fontpath,23);small=ImageFont.truetype(fontpath,17)
    log=(out/'ffmpeg.log').open('wb')
    encoder=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24','-video_size','960x600',
        '-framerate','30','-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','20','-pix_fmt','yuv420p',str(out/'comparison.mp4')],stdin=subprocess.PIPE,stdout=log,stderr=log)
    try:
        for frame in range(241):
            t=min(frame/30,8.);canvas=Image.new('RGB',(960,600),'#edf2f8');draw=ImageDraw.Draw(canvas)
            for j,(r,body) in enumerate(zip(chosen,bodies)):
                state=min(r['records'],key=lambda x:abs(x['time_s']-t))
                # Render recorded physical poses, never a prescribed spin angle.
                for hidden in bodies:p.resetBasePositionAndOrientation(hidden,[10,10,10],[0,0,0,1],physicsClientId=client)
                p.resetBasePositionAndOrientation(body,state['center_position_m'],state['principal_quaternion'],physicsClientId=client)
                view=p.computeViewMatrix([.095,-.12,.09],[0,0,.017],[0,0,1])
                proj=p.computeProjectionMatrixFOV(37,1,0.001,3)
                image=p.getCameraImage(480,480,viewMatrix=view,projectionMatrix=proj,renderer=p.ER_TINY_RENDERER,
                    lightDirection=[-1,-2,4],shadow=1,physicsClientId=client)
                canvas.paste(Image.fromarray(np.asarray(image[2],dtype=np.uint8).reshape(480,480,4)[:,:,:3]),(j*480,62))
                draw.text((j*480+18,14),('优化前：实心 Heart' if j==0 else '二次曲面优化后'),font=font,fill='#142944')
                draw.text((j*480+18,43),f'{props[r["scene"]["material"]][0]*1000:.2f} g  |  t = {t:.2f} s',font=small,fill='#425975')
                draw.text((j*480+18,548),f'倾角 {state["tilt_degrees"]:.2f}°   轴向转速 {state["axial_spin_hz"]:.2f} Hz',font=small,fill='#142944')
            draw.text((18,577),'PyBullet 实际轨迹 | 35 Hz 起转 · 2°扰动 · dt=0.00025 s | 外形相同、惯量不同',font=small,fill='#425975')
            if frame in (0,60,240):canvas.save(out/f'physics-{frame//30:02d}s.png')
            encoder.stdin.write(np.asarray(canvas).tobytes())
        encoder.stdin.close()
        if encoder.wait(timeout=45)!=0:raise RuntimeError('VIDEO_ENCODING_FAILED')
    finally:
        log.close();p.disconnect(client)

def finish(out,props,results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(12,4))
    for material,color in [('original','#d67538'),('optimized','#098f8b')]:
        for dt,style in [(.0005,'--'),(.00025,'-')]:
            r=next(r for r in results if r['scene']['material']==material and r['scene']['azimuth_degrees']==0 and r['scene']['dt_seconds']==dt)
            if r['status']=='FAILED':continue
            for ax,key in zip(axs,('tilt_degrees','axial_spin_hz')):
                ax.plot([x['time_s'] for x in r['records']],[x[key] for x in r['records']],style,color=color,label=f'{material}, dt={dt}')
    for ax,title in zip(axs,('0° direction: tilt / degrees','0° direction: axial spin / Hz')):
        ax.set(title=title,xlabel='Physical time / s');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'physics-curves.png',dpi=150);plt.close(fig)
    rows=''.join(f'<tr><td>{r["scene"]["material"]}</td><td>{r["scene"]["azimuth_degrees"]}°</td><td>{r["scene"]["dt_seconds"]}</td><td>{r.get("maximum_tilt_degrees",float("nan")):.3f}°</td><td>{r.get("final_spin_hz",float("nan")):.3f}</td><td>{r["status"]}</td></tr>' for r in results)
    html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>二次曲面 Heart · 物理对比</title>
<style>body{margin:0;background:#edf2f8;color:#142944;font:17px/1.65 system-ui,sans-serif}main{max-width:1180px;margin:auto;padding:35px 24px}h1{font-size:36px;line-height:1.2}h2{margin-top:40px}img,video{display:block;width:100%;border-radius:12px;margin:16px 0;background:white}.tag{display:inline-block;border-radius:20px;background:#d2eee9;padding:5px 16px}.note{background:white;border-left:4px solid #da8150;padding:14px 20px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:8px;border-bottom:1px solid #cbd5e1;text-align:left}.scroll{overflow:auto}a{color:#067e80}</style><main>
<span class="tag">实际优化结果 + 实际 PyBullet 轨迹</span><h1>外形不变，改变内部材料，观察旋转</h1>
<p>复用已有重建 Heart 的成功结果：<b>Faster 二次曲面与球面 Newton 求解</b>，目标为 <b>Spin-It 2014 惯量＋平方质心高度</b>。这是固定外形功能的二次曲面重实现；没有重新优化或换用体素结果。</p>
<h2>先看物理旋转</h2><video controls preload="metadata" poster="physics-00s.png" src="comparison.mp4"></video>
<p>视频方向提前固定为0°、时间步0.00025秒。画面播放求解器记录的物理姿态。起转35 Hz、2°倾斜、重力9.81、摩擦0.8，无持续驱动、无固定支点、阻尼及滚动／旋转摩擦为0。</p>
<div class="note">共同碰撞代理为原外形的凸包，已核对平面支撑距离；接触面分布与真实摩擦仍是近似。内腔浮点网格只展示，不作碰撞输入。质量、质心和完整非对角惯量来自保存的精确材料重新积分。此视频不是实物停止时间预测，也不是任意模型稳定保证。</div>
<h2>与 Spin-It 体素方案有什么区别？</h2><img src="method-comparison.png">
<div class="scroll"><table><tr><th></th><th>Spin-It 体素表示</th><th>本例 Faster 二次曲面表示</th></tr><tr><td>优化变量</td><td>各体素／八叉树叶子的填充率</td><td>一个全局二次多项式的10个归一化系数</td></tr><tr><td>内部边界</td><td>各单元共同确定材料</td><td>二次曲面与实际设计域求交</td></tr><tr><td>优化方法</td><td>高维 SLQP、自适应细分和二值化</td><td>球面 Newton、几何 Boolean、质量矩与导数</td></tr><tr><td>本页对照范围</td><td>图中方格是同一最终曲面的采样示意</td><td>实际历史优化材料与真实残差轨迹</td></tr></table></div>
<p><b>没有同输入的体素优化对照，因此本页不声称 Faster 比体素转得更久或计算更快。</b>相同最终惯量和外形会产生相同的刚体动力学；求解方法差异体现在如何得到内部材料。</p>
<h2>时间步检查与全部场景</h2><img src="physics-curves.png"><p>同一外形凸包，优化前后 × 八个扰动方向 × 两种时间步，共32个场景。全部保留，不挑最佳方向。没有完成第二种碰撞代理的独立验证。</p>
<div class="scroll"><table><tr><th>材料</th><th>方向</th><th>时间步</th><th>最大倾角</th><th>末态轴向Hz</th><th>状态</th></tr>'''+rows+'''</table></div>
<p><a href="material-properties.json">重新积分的材料属性与来源</a> · <a href="simulation-summary.json">全部场景结果</a> · <a href="simulation-protocol.json">固定仿真协议</a> · <a href="comparison.mp4">下载视频</a></p></main></html>'''
    (out/'index.html').write_text(html)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);a=parser.parse_args()
    out=a.output;out.mkdir(parents=True,exist_ok=False)
    began=time.monotonic();meshes,props,run,hashes=prepare(out)
    figures(out,meshes,props,run);results=simulate(out,props)
    render_video(out,props,results);finish(out,props,results)
    for path,sha in hashes.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=sha:raise ValueError('SOURCE_MODIFIED_DURING_DEMO')
    write(out/'complete.json',dict(status='VISUALIZATION_AND_CONTROLLED_PHYSICS_COMPLETE',elapsed_seconds=time.monotonic()-began,
        optimizer_invocations=0,scenarios=len(results),failed_scenarios=sum(r['status']=='FAILED' for r in results),
        paper_original=False,voxel_solver_comparison=False,manufacturing_approved=False,
        output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}))
    print('Completed',out,flush=True)

if __name__=='__main__':main()
