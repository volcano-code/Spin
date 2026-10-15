"""Paper-style views of saved qualified material; display transforms only."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

GRAY='#a7adb4'
PINK='#da71a6'
YELLOW='#e8bb40'
BLUE='#356fb0'


def circle_in_mm(record,frozen):
    scale=float(frozen['meters_per_model_unit'])*1000
    if not np.isfinite(scale) or scale<=0:raise ValueError('INVALID_FROZEN_LENGTH_UNIT')
    return np.asarray(record['circle_center'])*scale,record['circle_radius']*scale


def load_material(material, record):
    import trimesh
    from flower_render_physics import validate_material_meshes
    validate_material_meshes(material,record)
    meshes={}
    for name,info in record['files'].items():
        mesh=trimesh.load(material/info['file'],force='mesh',process=False)
        meshes[name]=(np.asarray(mesh.vertices),np.asarray(mesh.faces))
    return meshes


def scene_parts(record, mode):
    """Every part is the original saved mesh plus an explicit display offset."""
    fixed=list(record['rotor_layers'])
    if mode=='input':
        return [('design',PINK,1.,0.)]+[(n,GRAY,.16,0.) for n in fixed]
    if mode=='cutaway':
        return [('retained',YELLOW,1.,0.)]+[(n,GRAY,.7,0.) for n in fixed if n not in ('top','bottom')]
    if mode=='exploded':
        return [('retained',YELLOW,1.,0.),('top',GRAY,.7,.014),('bottom',GRAY,.7,-.014)]+[
            (n,GRAY,.6,0.) for n in fixed if n not in ('top','bottom')]+[
            (n,BLUE if n=='shaft' else GRAY,.8,.025 if n=='upper_cap' else -.025 if n=='lower_cap' else 0.)
            for n in record['clamp_layers']]
    if mode=='assembly':
        return [('retained',YELLOW,1.,0.)]+[(n,GRAY,.8,0.) for n in fixed]+[
            (n,GRAY,1.,0.) for n in record['clamp_layers']]
    raise ValueError('UNKNOWN_DISPLAY_MODE')


def render(folder,saved_check,physics,out,*,source_video=None,delivery_directory=None):
    from flower_result_figure import validate_status_sources
    status=validate_status_sources(folder,saved_material_check=saved_check,physics_summary=physics/'summary.json')
    material=Path(status['material_preparation_path']).parent
    record=json.loads((material/'output-preparation.json').read_text())
    frozen=json.loads((folder/'input.json').read_text())
    meshes=load_material(material,record)
    out.mkdir(parents=True,exist_ok=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.collections import PolyCollection
    font=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    def draw3(ax,mode):
        for name,color,opacity,dz in scene_parts(record,mode):
            v,f=meshes[name];points=v.copy();points[:,2]+=dz
            ax.add_collection3d(Poly3DCollection(points[f],facecolor=color,edgecolor='none',alpha=opacity))
        ax.set(xlim=(-.026,.026),ylim=(-.026,.026),zlim=(-.03,.032) if mode=='exploded' else (-.005,.007))
        ax.set_box_aspect((1,1,.8 if mode=='exploded' else .25),zoom=1.3)
        ax.view_init(elev=78 if mode=='cutaway' else 55,azim=-90);ax.set_axis_off()
    fig=plt.figure(figsize=(15,6),facecolor='white')
    titles=['固定层透明显示 / 粉色设计域','外层剖开 / 黄色为真实保留材料','完整三维装配']
    for k,mode in enumerate(('input','cutaway','assembly')):
        ax=fig.add_subplot(1,3,k+1,projection='3d');draw3(ax,mode)
        ax.set_title(titles[k],fontproperties=font,fontsize=15,pad=2)
    fig.text(.5,.95,'Spin-It Faster 花盘：输入 → 优化剖视 → 三维装配',ha='center',fontproperties=font,fontsize=20)
    fig.text(.5,.07,'灰：固定材料   粉：设计域   黄：优化保留材料   透明缺口：实际切除',ha='center',fontproperties=font,fontsize=12)
    fig.text(.5,.025,'重建方法验证；剖视外层隐藏仅用于展示。采用原优化与合格 OBJ，无新增求解。',ha='center',fontproperties=font,fontsize=10)
    fig.subplots_adjust(left=0,right=1,top=.84,bottom=.12,wspace=0)
    fig.savefig(out/'paper-three-columns.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,7),facecolor='white')
    front_parts=[('retained',YELLOW,1)]+[(name,GRAY,1) for name in record['rotor_layers'] if name not in ('top','bottom')]
    for name,color,alpha in front_parts:
        v,f=meshes[name];tri=v[f];top=tri[np.ptp(tri[:,:,2],axis=1)<1e-15]
        ax.add_collection(PolyCollection(top[:,:,:2]*1000,facecolor=color,edgecolor='none',alpha=alpha))
    from matplotlib.patches import Circle
    center,radius=circle_in_mm(record,frozen)
    ax.add_patch(Circle(center,radius,fill=False,color='#555',ls='--',lw=1.3))
    ax.plot(0,0,'+',color=BLUE,ms=12);ax.annotate('实际右下侧切除',xy=center,xytext=(12,-22),
        fontproperties=font,arrowprops={'arrowstyle':'->','color':'#555'})
    ax.set(xlim=(-27,27),ylim=(-27,27),xlabel='+x → [mm]',ylabel='+y ↑ [mm]');ax.set_aspect('equal')
    ax.set_title('中间层正面：黄色保留，灰色固定，白色为空腔或原始镂空',fontproperties=font,fontsize=13)
    fig.tight_layout();fig.savefig(out/'fixed-front.png',dpi=160);plt.close(fig)
    fig=plt.figure(figsize=(8,8),facecolor='white');ax=fig.add_subplot(projection='3d');draw3(ax,'exploded')
    ax.set_title('三层展开与独立夹持组件',fontproperties=font,fontsize=18)
    fig.text(.5,.04,'仅展示变换；仿真始终使用原装配位置。',fontproperties=font,ha='center')
    fig.savefig(out/'exploded.png',dpi=150);plt.close(fig)
    import plotly.graph_objects as go
    traces=[]
    modes=('input','cutaway','assembly','exploded')
    counts=[]
    for mode in modes:
        parts=scene_parts(record,mode);counts.append(len(parts))
        for name,color,opacity,dz in parts:
            v,f=meshes[name]
            traces.append(go.Mesh3d(x=v[:,0]*1000,y=v[:,1]*1000,z=(v[:,2]+dz)*1000,
                i=f[:,0],j=f[:,1],k=f[:,2],color=color,opacity=opacity,name=name,
                visible=mode=='cutaway',flatshading=True,hoverinfo='name'))
    starts=np.cumsum([0]+counts)
    labels=['输入','优化剖视','完整装配','三层展开']
    buttons=[dict(label=label,method='update',args=[{'visible':[starts[k]<=i<starts[k+1] for i in range(len(traces))]}])
             for k,label in enumerate(labels)]
    plot=go.Figure(traces);plot.update_layout(height=670,margin=dict(l=0,r=0,t=50,b=0),paper_bgcolor='white',
        scene=dict(aspectmode='data',xaxis_title='x [mm]',yaxis_title='y [mm]',zaxis_title='z [mm]',
                   camera=dict(eye=dict(x=0,y=-1.4,z=2.2),up=dict(x=0,y=0,z=1))),
        updatemenus=[dict(type='buttons',direction='right',x=0,y=1.1,buttons=buttons)])
    fragment=plot.to_html(full_html=False,include_plotlyjs=True)
    video=(Path(source_video) if source_video else physics.parent.parent/'render-job-01/output/flower-physical-comparison.mp4').resolve()
    video_manifest=video.parent/'render-manifest.json'
    motion_identity=json.loads(video_manifest.read_text())
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    if (motion_identity['source_physics_sha256']!=sha(physics/'summary.json') or
        motion_identity['video_sha256']!=sha(video) or motion_identity['selected_phase']!=0 or
        motion_identity['selected_dt_seconds']!=.00025):raise ValueError('SOURCE_VIDEO_PHYSICS_IDENTITY_MISMATCH')
    import os
    video_link='paper-physical-comparison.mp4'
    if delivery_directory is not None:
        delivery=Path(delivery_directory)
        package=json.loads((delivery/'delivery-manifest.json').read_text())
        if package['material_preparation_sha256']!=sha(material/'output-preparation.json') or package['saved_check_sha256']!=sha(saved_check):
            raise ValueError('DELIVERY_MATERIAL_IDENTITY_MISMATCH')
        for name,digest in package['archives'].items():
            if sha(delivery/name)!=digest:raise ValueError('DELIVERY_ARCHIVE_CHANGED')
        zip_link=os.path.relpath(delivery/'final-assembly-obj.zip',out)
        diagnostic_link=os.path.relpath(delivery/'diagnostic-obj.zip',out)
        downloads='<a href="'+zip_link+'" download>下载最终装配 OBJ 与物理属性</a><p><a href="'+diagnostic_link+'" download>下载原始／设计域／空腔诊断包</a></p>'
    else:
        legacy=folder/'qualified-flower-obj-v1.zip'
        if not legacy.exists():raise ValueError('QUALIFIED_DELIVERY_PACKAGE_REQUIRED')
        downloads='<a href="'+os.path.relpath(legacy,out)+'" download>下载历史材料诊断包</a>'
    html='''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
    <title>花盘真实优化与三层结构</title><style>body{max-width:1280px;margin:32px auto;padding:0 24px;color:#22313d;font:16px system-ui;background:#fff}h1{font-size:28px}img,video{width:100%;height:auto}section{margin:34px 0}a{color:#25649b}.pair{display:grid;grid-template-columns:1fr 1fr;gap:24px}.note{color:#64717b}button{cursor:pointer}</style>
    <h1>Spin-It Faster：多瓣花盘真实优化</h1><p>灰色固定层、粉色设计域、黄色保留材料。右下侧切除来自实际 disk / dogleg 求解。</p>
    <p class="note">RECONSTRUCTED_METHOD_VALIDATION · 合格 OBJ · 64 个实际轴承／柔性支撑场景。48个原有孔保持；仅中间层新增圆形切除。</p>
    <section><img src="paper-three-columns.png" alt="论文式三栏实际材料图"></section>
    <section class="pair"><img src="fixed-front.png" alt="固定正面右下侧切除"><img src="exploded.png" alt="三层展开"></section>
    <section><h2>交互装配</h2><p>拖动旋转，选择输入、剖视或展开；坐标单位为毫米。</p>'''+fragment+'''</section>
    <section><h2>真实物理对比</h2><p>先结合上方固定正面识别缺口，再观察旋转。视频读取保存轨迹，采用首个登记相位；含正常速度及0.1倍慢放。</p>
    <video controls preload="metadata" src="'''+video_link+'''"></video><p class="note">理想轴承与未标定柔性支撑工程模型，优化前后共用支撑参数。</p></section>
    <section>'''+downloads+'''</section></html>'''
    (out/'index.html').write_text(html,encoding='utf-8')
    manifest=dict(schema='qualified-flower-paper-display-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',
        **status,circle_center=record['circle_center'],circle_radius=record['circle_radius'],
        coordinate_contract='MODEL_X_RIGHT_Y_UP_FIXED_FRONT; meshes_and_axes_in_mm',
        mesh_sha256={name:info['sha256'] for name,info in record['files'].items()},
        source_video_manifest_path=str(video_manifest),source_video_manifest_sha256=sha(video_manifest),
        meters_per_model_unit=frozen['meters_per_model_unit'],
        display_transforms={mode:scene_parts(record,mode) for mode in modes},
        optimizer_invocations=0,physical_simulator_invocations=0,
        delivery_manifest_sha256=sha(Path(delivery_directory)/'delivery-manifest.json') if delivery_directory else None,
        renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'display-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Preserve the recorded motion frames; add a static front-view sidebar and
    # a three-second explanatory introduction without new integration.
    from PIL import Image,ImageDraw,ImageFont
    panel=Image.new('RGB',(320,800),'white');draw=ImageDraw.Draw(panel)
    small=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',18)
    front=Image.open(out/'fixed-front.png').convert('RGB');front.thumbnail((320,400))
    panel.paste(front,(0,55))
    for i,line in enumerate(['固定模型正面','右下侧为实际切除','运动画面读取真实轨迹','灰色：固定层','黄色：保留材料','正常速度＋0.1倍慢放','轴承／柔性支撑模型']):
        draw.text((15,12 if i==0 else 455+32*(i-1)),line,font=small,fill='#263746')
    panel.save(out/'video-sidebar.png')
    intro=Image.new('RGB',(1280,800),'white');main=Image.open(out/'paper-three-columns.png').convert('RGB')
    main.thumbnail((1280,720));intro.paste(main,((1280-main.width)//2,(800-main.height)//2));intro.save(out/'video-intro.png')
    import subprocess
    command=['ffmpeg','-nostdin','-n','-threads','1','-loop','1','-framerate','30','-t','3','-i',str(out/'video-intro.png'),
        '-i',str(video),'-loop','1','-framerate','30','-i',str(out/'video-sidebar.png'),
        '-filter_complex','[1:v]pad=1280:800:0:0:white[base];[base][2:v]overlay=960:0:shortest=1[motion];[0:v]setsar=1[intro];[intro][motion]concat=n=2:v=1:a=0[out]',
        '-map','[out]','-an','-c:v','libx264','-threads','1','-preset','fast','-crf','20','-pix_fmt','yuv420p',str(out/video_link)]
    with (out/'video.log').open('xb') as log:
        subprocess.run(command,stdout=log,stderr=log,check=True,timeout=600)
    manifest.update(source_video_sha256=hashlib.sha256(video.read_bytes()).hexdigest(),
        video_sha256=hashlib.sha256((out/video_link).read_bytes()).hexdigest(),
        intro_display_seconds=3,static_sidebar=True)
    (out/'display-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('folder','saved-check','physics','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--source-video',type=Path)
    parser.add_argument('--delivery-directory',type=Path)
    args=parser.parse_args();render(args.folder,args.saved_check,args.physics,args.output,source_video=args.source_video,delivery_directory=args.delivery_directory)
