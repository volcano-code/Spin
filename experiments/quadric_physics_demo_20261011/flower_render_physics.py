"""Render saved admitted physical trajectories; never integrate or invent motion."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np


def validate_material_meshes(material,record):
    for info in record['files'].values():
        path=material/info['file']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=info['sha256']:
            raise ValueError('RENDER_MATERIAL_MESH_IDENTITY_MISMATCH:'+info['file'])


def render(material, physics, out, *, paper_colors=False):
    import pybullet as p
    from PIL import Image, ImageDraw, ImageFont
    record=json.loads((material/'output-preparation.json').read_text())
    summary=json.loads((physics/'summary.json').read_text())
    if summary['recorded_scenes'] != summary['expected_scenes'] or summary['expected_scenes'] != 64:
        raise ValueError('INCOMPLETE_PHYSICS_MATRIX')
    if summary['sources']['material_preparation_sha256'] != hashlib.sha256((material/'output-preparation.json').read_bytes()).hexdigest():
        raise ValueError('PHYSICS_MATERIAL_IDENTITY_MISMATCH')
    validate_material_meshes(material,record)
    out.mkdir(parents=True,exist_ok=False)
    trajectories={}
    for model in ('ideal','flexible'):
        for variant in ('original','optimized'):
            name=f'{model}-{variant}-phase0-dt0.00025'
            entry=next(x for x in summary['scenes'] if x['id']==name)
            path=physics/(name+'.npz')
            if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['trajectory_sha256']:
                raise ValueError('PHYSICS_TRAJECTORY_HASH_MISMATCH')
            trajectories[(model,variant)]=np.load(path)
    client=p.connect(p.DIRECT)
    try:
        bodies={}
        clamp_names=list(record['clamp_layers'])
        for variant in ('original','optimized'):
            names=['design' if variant=='original' else 'retained',*record['rotor_layers'],*clamp_names]
            bodies[variant]=[]
            for name in names:
                path=material/(name+'.obj')
                is_clamp=name in clamp_names
                color=([.16,.37,.74,1] if is_clamp else [.60,.65,.69,1] if name in ('bottom','top') else
                       [.82,.60,.17,1] if name=='outer-race' else [.85,.32,.57,1] if variant=='original' else [.05,.68,.52,1])
                if paper_colors:
                    color=([.64,.68,.72,1] if is_clamp or name in record['rotor_layers'] else
                           [.85,.44,.65,1] if variant=='original' else [.91,.73,.25,1])
                shape=p.createVisualShape(p.GEOM_MESH,fileName=str(path),rgbaColor=color,physicsClientId=client)
                body=p.createMultiBody(baseMass=0,baseVisualShapeIndex=shape,baseCollisionShapeIndex=-1,
                                       basePosition=[10,0,0],physicsClientId=client)
                bodies[variant].append((body,is_clamp))
        width,height=480,340
        view=p.computeViewMatrix([.032,-.060,.059],[0,0,0],[0,0,1])
        projection=p.computeProjectionMatrixFOV(43,width/height,.001,2)
        fontpath='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
        font=ImageFont.truetype(fontpath,19);small=ImageFont.truetype(fontpath,16)
        with (out/'video.log').open('xb') as log:
            process=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24',
                '-video_size','960x800','-framerate','30','-i','pipe:0','-an','-c:v','libx264',
                '-preset','fast','-crf','20','-pix_fmt','yuv420p',str(out/'flower-physical-comparison.mp4')],
                stdin=subprocess.PIPE,stdout=log,stderr=log)
            try:
                for frame in range(540):
                    slow=frame>=240;t=(frame-240)/300 if slow else frame/30
                    image=Image.new('RGB',(960,800),'#f1f4f7');draw=ImageDraw.Draw(image)
                    for row,model in enumerate(('ideal','flexible')):
                        for col,variant in enumerate(('original','optimized')):
                            data=trajectories[(model,variant)];index=min(int(round(t/.00025)),len(data['t'])-1)
                            angle=float(data['angle'][index]);q=p.getQuaternionFromEuler([0,0,angle])
                            xy=data['support_xy_m'][index];position=[float(xy[0]),float(xy[1]),0]
                            for items in bodies.values():
                                for body,_ in items:p.resetBasePositionAndOrientation(body,[10,0,0],[0,0,0,1],physicsClientId=client)
                            for body,is_clamp in bodies[variant]:
                                p.resetBasePositionAndOrientation(body,position,[0,0,0,1] if is_clamp else q,physicsClientId=client)
                            pixels=p.getCameraImage(width,height,view,projection,renderer=p.ER_TINY_RENDERER,
                                lightDirection=[-2,-3,5],physicsClientId=client)[2]
                            image.paste(Image.fromarray(np.asarray(pixels,np.uint8).reshape(height,width,4)[:,:,:3]),(col*480,40+row*375))
                            label=('理想轴承' if model=='ideal' else '柔性支撑')+' / '+('优化前' if variant=='original' else '优化后')
                            draw.text((15+col*480,10+row*375),label,font=font,fill='#15304b')
                            draw.text((15+col*480,354+row*375),
                                f'转速 {data["omega"][index]/(2*np.pi):.2f} Hz · 支撑位移 {np.linalg.norm(xy)*1e6:.2f} μm',font=small,fill='#15304b')
                    draw.text((15,758),f'实际物理轨迹 t={t:.3f}s / '+('0.1倍慢放' if slow else '正常速度；高速运动有帧采样效应'),font=small,fill='#15304b')
                    caption=('灰色固定件，粉色原始设计材料，黄色优化材料。首相位；未标定工程支撑。' if paper_colors else
                             '蓝色为夹持件；首个登记相位，无马达驱动。柔性支撑为未实测标定工程模型。')
                    draw.text((15,780),caption,font=small,fill='#15304b')
                    if frame==0:image.save(out/'preview.png')
                    process.stdin.write(np.asarray(image).tobytes())
            finally:
                process.stdin.close()
            if process.wait(timeout=60):raise RuntimeError('VIDEO_ENCODING_FAILED')
        (out/'render-manifest.json').write_text(json.dumps(dict(schema='saved-physics-render-v1',
            source_physics_sha256=hashlib.sha256((physics/'summary.json').read_bytes()).hexdigest(),
            video_sha256=hashlib.sha256((out/'flower-physical-comparison.mp4').read_bytes()).hexdigest(),
            selected_phase=0,selected_dt_seconds=.00025,physical_simulator_invocations=0,
            normal_physical_interval=[0,8],slow_physical_interval=[0,1],slow_playback_ratio=.1,
            paper_colors=paper_colors),indent=2)+'\n')
    finally:p.disconnect(client)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('material','physics','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--paper-colors',action='store_true')
    args=parser.parse_args();render(args.material,args.physics,args.output,paper_colors=args.paper_colors)
