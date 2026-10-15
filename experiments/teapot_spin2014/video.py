"""Render stored dynamics, never prescribe an angular-speed animation."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.spatial.transform import Rotation, Slerp
from PIL import Image, ImageDraw, ImageFont
from teapot import STRATEGIES, write, file_hashes, check_files
from admission import check_video_context


def comparison_scope(names):
    names=set(names)
    if names=={'original'}:return 'ORIGINAL_BASELINE_ONLY'
    if names=={'original',*STRATEGIES}:return 'FOUR_WAY_RECONSTRUCTED_COMPARISON'
    return 'PARTIAL_COMPARISON'


def trajectory_pose(record, t, center):
    if 'samples' in record:
        times=np.array([s['time'] for s in record['samples']])
        qs=np.array([s['model_quaternion'] for s in record['samples']])
        positions=np.array([s['position'] for s in record['samples']])
    else:
        times=np.asarray(record['time']);qs=np.asarray(record['quaternions']);positions=None
    if (not np.isfinite(times).all() or not np.isfinite(qs).all() or not np.isfinite(t)
        or len(times)<2 or np.any(np.diff(times)<=0) or t<times[0] or t>times[-1]):
        raise ValueError('invalid or out-of-range trajectory')
    R=Slerp(times,Rotation.from_quat(qs))(t)
    if positions is None:return np.zeros(3),R
    if not np.isfinite(positions).all():raise ValueError('nonfinite positions')
    center=np.asarray(center,float)
    if center.shape!=(3,) or not np.isfinite(center).all():raise ValueError('invalid center')
    position=np.array([np.interp(t,times,positions[:,i]) for i in range(3)])-R.apply(center)
    return position,R


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--physics',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--hz',type=float,default=35.);p.add_argument('--profile',default='default')
    a=p.parse_args();a.output.mkdir(exist_ok=False)
    manifest=json.loads((a.input/'manifest.json').read_text());check_files(manifest['geometry_hashes'])
    import pybullet as b
    client=b.connect(b.DIRECT)
    try:
        shape=b.createVisualShape(b.GEOM_MESH,fileName=str(a.input/'outer.obj'),meshScale=[.025]*3,
            rgbaColor=[.3,.63,.46,1],physicsClientId=client)
        body=b.createMultiBody(baseMass=0,baseVisualShapeIndex=shape,physicsClientId=client)
        view=b.computeViewMatrix([.07,-.12,.035],[0,0,.021],[0,0,1])
        projection=b.computeProjectionMatrixFOV(38,4/3,.001,2.)
        def snap():
            im=b.getCameraImage(512,384,viewMatrix=view,projectionMatrix=projection,
                lightDirection=[-2,-3,5],shadow=1,renderer=b.ER_TINY_RENDERER,physicsClientId=client)
            return Image.fromarray(np.asarray(im[2],dtype=np.uint8).reshape(384,512,4)[:,:,:3])
        snap().save(a.output/'teapot-input.png')
        if a.physics is None:
            write(a.output/'scope.json',dict(scope='INPUT_GEOMETRY_PREVIEW_NOT_SIMULATION',source=file_hashes([a.input/'outer.obj'])))
            return
        report=json.loads((a.physics/'report.json').read_text());check_files(json.loads((a.physics/'source.json').read_text())['files'])
        protocol=report['protocol']; profile='ideal' if protocol=='ideal' else a.profile
        titles={'original':'Original solid','inertia':'Inertia only','lowering':'Inertia + lowering',
            'mass_lowering':'Inertia + mass/height'}
        names=['original',*STRATEGIES];records={};centers={m['label']:m['center_m'] for m in report['models']}
        paths=[a.physics/'report.json',a.physics/'source.json',a.input/'manifest.json',Path(__file__)]
        for name in names:
            path=a.physics/f'{name}-{a.hz:g}hz-{profile}-dir0.json'
            if path.is_file():
                check_video_context(report,manifest['problem_digest'],path)
                records[name]=json.loads(path.read_text());paths.append(path)
        if not records:raise ValueError('NO_ACTUAL_TRAJECTORIES')
        evidence=file_hashes(paths);evidence.update(json.loads((a.physics/'source.json').read_text())['files'])
        evidence.update(manifest['geometry_hashes']);evidence.update(report['trajectory_hashes'])
        evidence.update(file_hashes([Path(__file__).with_name('admission.py')]))
        write(a.output/'source.json',dict(files=evidence))
        if protocol=='contact':
            floor_shape=b.createVisualShape(b.GEOM_BOX,halfExtents=[.4,.4,.02],rgbaColor=[.70,.66,.58,1],physicsClientId=client)
            b.createMultiBody(baseMass=0,baseVisualShapeIndex=floor_shape,basePosition=[0,0,-.02],physicsClientId=client)
        font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',17)
        small=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',14)
        colors=[[.30,.52,.80,1],[.93,.65,.18,1],[.45,.68,.36,1],[.80,.38,.25,1]]
        with (a.output/'ffmpeg.log').open('x') as log:
            process=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24',
                '-video_size','1024x876','-framerate','30','-i','pipe:0','-an','-c:v','libx264',
                '-preset','fast','-crf','20','-pix_fmt','yuv420p',str(a.output/'teapot-comparison.mp4')],
                stdin=subprocess.PIPE,stdout=log,stderr=log)
            try:
                for frame in range(241):
                    t=frame/30; canvas=Image.new('RGB',(1024,876),'#f4f3ef');draw=ImageDraw.Draw(canvas)
                    mode='IDEAL FIXED PIVOT (no table collision)' if protocol=='ideal' else 'PYBULLET FREE TABLE CONTACT'
                    draw.text((20,12),f'Utah Teapot | {mode}',font=font,fill='#20252b')
                    scope=comparison_scope(records).replace('_',' ')
                    draw.text((20,38),f'{scope} | {a.hz:g} rps | t={t:.2f}s',font=small,fill='#40454b')
                    for i,name in enumerate(names):
                        x=(i%2)*512;y=66+(i//2)*384
                        if name in records:
                            origin,R=trajectory_pose(records[name],t,centers[name])
                            b.resetBasePositionAndOrientation(body,origin.tolist(),R.as_quat().tolist(),physicsClientId=client)
                            b.changeVisualShape(body,-1,rgbaColor=colors[i],physicsClientId=client)
                            canvas.paste(snap(),(x,y))
                            tilt=float(np.rad2deg(np.arccos(np.clip(R.as_matrix()[2,2],-1,1))))
                            text=f'{titles[name]} | tilt {tilt:.2f} deg'
                        else:text=f'{titles[name]} | UNSEALED / NO TRAJECTORY'
                        draw.rectangle((x,y,x+511,y+27),fill='#e7e5de');draw.text((x+12,y+4),text,font=small,fill='#18202b')
                    draw.text((20,843),'Actual stored trajectories, 1x time. Numerical observations; not author-original or physical validation.',font=small,fill='#40454b')
                    if frame in [0,120,240]:canvas.save(a.output/f'frame-{frame:03d}.png')
                    process.stdin.write(np.asarray(canvas).tobytes())
                    if frame%60==0:print('VIDEO FRAME',frame,'/240',flush=True)
            finally:process.stdin.close()
            if process.wait(timeout=60)!=0:raise RuntimeError('video encoding failed')
        check_files(evidence)
        write(a.output/'video.json',dict(frames=241,fps=30,simulation_start=0,simulation_end=8,playback_speed=1,
            display='SOFTWARE_RENDERED_ACTUAL_TRAJECTORY_REPLAY_NOT_LIVE_DESKTOP',protocol=protocol,
            selected_direction=0,selected_hz=a.hz,profile=profile,
            all_directions_in_report=all({r.get('direction') for r in report['runs'] if r.get('label')==name
                and r.get('hz')==a.hz and r.get('profile')==profile and r['status']=='RECORDED'}==set(range(8)) for name in records),
            scope='RECONSTRUCTED_NOT_ORIGINAL_PAPER',comparison_scope=comparison_scope(records),
            files=file_hashes([a.output/'teapot-comparison.mp4'])))
    finally:b.disconnect(client)


if __name__=='__main__':main()
