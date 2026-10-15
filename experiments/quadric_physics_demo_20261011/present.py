"""Present all time-step outcomes of saved physics; no new dynamics or fitting."""
import hashlib
import json
import subprocess
from pathlib import Path
import sys
import numpy as np
from PIL import Image,ImageDraw,ImageFont

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DATA=HERE/'job-v1/output'
OUT=HERE/'presentation-v1'
sys.path.insert(0,str(ROOT/'experiments/heart_stem_contact'))

def main():
    import pybullet as p
    from bullet_contact_common_frame import create_body
    OUT.mkdir(exist_ok=False)
    props_data=json.loads((DATA/'material-properties.json').read_text())
    props={name:(v['mass_kg'],np.array(v['center_m']),np.array(v['inertia_kg_m2']),
                 np.array(v['principal_moments']),np.array(v['principal_axes'])) for name,v in props_data['materials'].items()}
    rows=[];inputs={}
    for dt in (.0005,.00025):
        for name in ('original','optimized'):
            path=DATA/f'{name}-default-dir000-dt{dt:g}.json'
            inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
            row=json.loads(path.read_text());assert row['status']!='FAILED';rows.append(row)
    client=p.connect(p.DIRECT)
    floor=p.createVisualShape(p.GEOM_BOX,halfExtents=[.4,.4,.001],rgbaColor=[.88,.91,.96,1],physicsClientId=client)
    p.createMultiBody(baseMass=0,baseVisualShapeIndex=floor,basePosition=[0,0,-.001],physicsClientId=client)
    bodies=[]
    for name,color in [('original',(.91,.48,.25,1)),('optimized',(.02,.62,.6,1))]:
        bodies.append(create_body(p,client,DATA/'outer-si.obj',[DATA/'collision-hull-si.obj'],props[name],np.eye(3),np.zeros(3),35,color))
    fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
    title=ImageFont.truetype(fp,24);font=ImageFont.truetype(fp,19);small=ImageFont.truetype(fp,16)
    with (OUT/'ffmpeg.log').open('wb') as log:
        encoder=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24','-video_size','960x840','-framerate','30','-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','20','-pix_fmt','yuv420p',str(OUT/'comparison-all-timesteps.mp4')],stdin=subprocess.PIPE,stdout=log,stderr=log)
        for frame in range(241):
            t=min(frame/30,8);canvas=Image.new('RGB',(960,840),'#edf2f8');draw=ImageDraw.Draw(canvas)
            draw.text((16,6),'相同外形、不同材料：PyBullet 实际旋转对比',font=title,fill='#142944')
            draw.text((16,36),'上下两种时间步均展示 · 接触结果尚未收敛，不能据此认定稳定 8 秒',font=small,fill='#a64425')
            for k,row in enumerate(rows):
                j=k%2;offset=60+(k//2)*370;state=min(row['records'],key=lambda x:abs(x['time_s']-t))
                for body in bodies:p.resetBasePositionAndOrientation(body,[10,10,10],[0,0,0,1],physicsClientId=client)
                p.resetBasePositionAndOrientation(bodies[j],state['center_position_m'],state['principal_quaternion'],physicsClientId=client)
                view=p.computeViewMatrix([.095,-.12,.09],[0,0,.017],[0,0,1]);projection=p.computeProjectionMatrixFOV(37,1.6,.001,3)
                image=p.getCameraImage(480,300,viewMatrix=view,projectionMatrix=projection,renderer=p.ER_TINY_RENDERER,lightDirection=[-1,-2,4],shadow=1,physicsClientId=client)
                canvas.paste(Image.fromarray(np.asarray(image[2],dtype=np.uint8).reshape(300,480,4)[:,:,:3]),(j*480,offset+40))
                draw.text((j*480+14,offset+4),('原始实心' if j==0 else '二次曲面优化')+f'   dt={row["scene"]["dt_seconds"]} s',font=font,fill='#142944')
                draw.text((j*480+14,offset+341),f'倾角 {state["tilt_degrees"]:.1f}°   轴向转速 {state["axial_spin_hz"]:.2f} Hz',font=small,fill='#142944')
            draw.text((16,808),f't = {t:.2f} s  |  35 Hz 起转 · 2°扰动 · 无驱动 · 共同外形凸包接触代理',font=font,fill='#425975')
            if frame in (0,60,240):canvas.save(OUT/f'comparison-{frame//30:02d}s.png')
            encoder.stdin.write(np.asarray(canvas).tobytes())
        encoder.stdin.close();assert encoder.wait(timeout=45)==0
    p.disconnect(client)
    cutaway()
    summary=json.loads((DATA/'simulation-summary.json').read_text())
    table=''.join(f'<tr><td>{r["scene"]["material"]}</td><td>{r["scene"]["azimuth_degrees"]}°</td><td>{r["scene"]["dt_seconds"]}</td><td>{r["maximum_tilt_degrees"]:.2f}°</td><td>{r["final_spin_hz"]:.3f} Hz</td></tr>' for r in summary['results'])
    html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Heart 二次曲面优化与物理旋转</title><style>body{margin:0;background:#edf2f8;color:#142944;font:17px/1.7 system-ui,sans-serif}main{max-width:1120px;margin:auto;padding:30px 24px}h1{font-size:36px;line-height:1.25}h2{margin-top:38px}img,video{display:block;width:100%;border-radius:12px;margin:14px 0}.notice{border-left:5px solid #cb6735;background:#fff2e5;padding:16px 22px}.pill{background:#ccece7;padding:6px 15px;border-radius:18px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border-bottom:1px solid #c1ccda;padding:9px;text-align:left}a{color:#067e80}</style><main><span class="pill">已有真实优化结果 · 新执行32个物理场景</span><h1>用 Faster 的方法，实现 Spin-It 的固定外形功能</h1><p>直接复用已成功的重建 Heart：二次曲面表示、CGAL Boolean、球面 Newton；目标采用 Spin-It 2014 的惯量＋平方质心高度。保持原外形，不增添握柄或改变材料。</p>
<h2>真实旋转：同时看两种时间步</h2><video controls preload="metadata" poster="comparison-00s.png" src="comparison-all-timesteps.mp4"></video><div class="notice"><b>时间步敏感，接触仿真尚未收敛。</b>0.00025秒时优化体在8秒末仍约33.8 Hz；0.0005秒时部分优化场景倾倒。这里公开两种结果，不能把较细时间步的画面当作稳定8秒的证明。视频来自真实轨迹，方向0°在运行前固定。</div><p>原始质量18.165 g，优化后15.450 g，减重14.94%。质心高度25.715 → 25.849 mm；这个复合目标的结果并没有让质心单独下降。仿真分别加载重新积分的质量、质心和完整非对角惯量。</p>
<h2>内部到底如何改变？</h2><img src="cutaway.png"><p>灰色为固定外壳；青色为实际优化保留材料；橙色为空腔。剖切仅供展示。权威材料使用保存的有理坐标，内部浮点网格没有作为碰撞体。</p>
<h2>这与 Spin-It 体素求解有什么差别？</h2><img src="../job-v1/output/method-comparison.png"><div class="scroll"><table><tr><th>比较</th><th>Spin-It 体素路径</th><th>本例二次曲面路径</th></tr><tr><td>变量</td><td>各体素／八叉树叶子的填充率</td><td>全局二次多项式：10个归一化系数</td></tr><tr><td>材料边界</td><td>单元填充共同决定</td><td>三角化二次曲面与设计域实际求交</td></tr><tr><td>求解</td><td>SLQP、自适应细分、二值化</td><td>球面 Newton、Boolean、质量矩与导数</td></tr><tr><td>本页证据</td><td>方格为同一曲面的采样示意，未运行体素优化</td><td>已有真实优化材料与残差轨迹</td></tr></table></div><p>平滑线为解析二次曲面截面，实际保留材料为离散 Boolean 结果，两者不被视作精确相等。没有同输入的体素优化对照，所以不作转动性能或加速比结论。刚体运动由几何与质量属性决定，不会因求解器名字不同而改变。</p>
<h2>完整场景与可复查数据</h2><img src="../job-v1/output/physics-curves.png"><p>优化前后 × 八方向 × 两时间步＝32个场景，全部保存。共同外形凸包保留了平面支撑位置，接触面分布仍是代理近似；尚未做第二碰撞代理对照。摩擦0.8，初始35 Hz、2°倾斜；重力9.81，无持续驱动、无固定支点、阻尼和滚动／旋转摩擦为0。没有实测标定，不预测现实停止时间。</p><div class="scroll"><table><tr><th>材料</th><th>扰动方向</th><th>dt</th><th>最大倾角</th><th>8秒轴向转速</th></tr>'''+table+'''</table></div><p><a href="comparison-all-timesteps.mp4">下载完整对比视频</a> · <a href="../job-v1/output/material-properties.json">材料属性及来源</a> · <a href="../job-v1/output/simulation-summary.json">32场景统计</a> · <a href="../job-v1/output/simulation-protocol.json">仿真协议</a></p></main></html>'''
    (OUT/'index.html').write_text(html)
    for path,digest in inputs.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
    (OUT/'manifest.json').write_text(json.dumps(dict(scope='PRESENTATION_OF_SAVED_PHYSICS',new_simulations=0,
        timestep_convergence='NOT_PASSED',source_trajectories=inputs,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file()}),indent=2))

def cutaway():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import importlib.util
    spec=importlib.util.spec_from_file_location('archived_render',ROOT/'experiments/faster_quadric_visualization_20261011/render.py')
    renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)
    folder=Path(json.loads((DATA/'material-properties.json').read_text())['source'])/'material-geometry'
    manifest=json.loads((folder/'manifest.json').read_text())
    meshes={k:(np.load(folder/v['vertices']['name']),np.load(folder/v['faces']['name'])) for k,v in manifest['meshes'].items()}
    fig=plt.figure(figsize=(14,5),facecolor='#f4f7fc')
    for j,(names,title) in enumerate([(['fixed','retained'],'固定壳＋保留材料：正面剖切'),(['retained'],'实际保留材料'),(['cavity'],'实际移除空腔')]):
        ax=fig.add_subplot(1,3,j+1,projection='3d',computed_zorder=False)
        for name in names:
            v,f=meshes[name];tris=renderer.clip_triangles(v,f,axis=1,value=0) if j==0 else v[f]
            color={'fixed':'#b6c4d5','retained':'#0c9a94','cavity':'#e89b4f'}[name]
            ax.add_collection3d(Poly3DCollection(tris,facecolors=color,edgecolors='none',
                alpha=.22 if name=='fixed' else 1,zorder=1 if name=='fixed' else 2))
        ax.set(xlim=(-1,1),ylim=(-.48,.48),zlim=(0,1.9),title=title)
        ax.set_box_aspect((2,.96,1.9));ax.view_init(elev=15,azim=73);ax.set_axis_off()
    fig.text(.04,.035,'实际保存的三角面，仅调整观察方向与展示剖切；不平滑、不修补、不参与碰撞计算。',fontsize=12)
    fig.subplots_adjust(left=.02,right=.98,top=.91,bottom=.09,wspace=.01);fig.savefig(OUT/'cutaway.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
