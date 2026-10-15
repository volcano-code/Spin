"""Serial approved flower engine/render work and noodles cost diagnosis."""
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import zipfile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job


def main():
    out=HERE/'paper-flower-method-v2/mujoco-blender-v1';out.mkdir(exist_ok=False)
    bundle=HERE/'paper-flower-method-v2/table-scene-v1'
    coordinator=1766941
    events=(out/'events.jsonl').open('x')
    def event(status,**data):
        row=dict(status=status,**data);events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno());print(json.dumps(row),flush=True)
    def signal_coordinator(sig):
        path=Path(f'/proc/{coordinator}/cmdline')
        if path.exists() and b'rounded_production_sequence.py' in path.read_bytes():os.kill(coordinator,sig)
        else:raise RuntimeError('HEART_COORDINATOR_IDENTITY_CHANGED')
    signal_coordinator(signal.SIGSTOP)
    try:
        # Read after stopping the dispatcher to include any point it just started.
        rows=[json.loads(x) for x in (ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-production-sequence-01/events.jsonl').read_text().splitlines()]
        current=Path(next(r['job'] for r in reversed(rows) if r['status']=='POINT_STARTED'))
        event('WAIT_ACTIVE_HEART_POINT',job=str(current))
        deadline=time.monotonic()+3660
        while not (current/'complete.json').exists():
            if time.monotonic()>deadline:raise RuntimeError('ACTIVE_HEART_JOB_NOT_FINISHED')
            time.sleep(2)
        python=ROOT/'.venv-flower-mujoco/bin/python'
        package=ROOT/'.venv-flower-mujoco/lib/python3.13/site-packages/mujoco'
        files={'worker':HERE/'flower_mujoco.py','support':HERE/'flexible_bearing.py','config':bundle/'table-scene.json'}
        files.update({str(p.relative_to(bundle)):p for p in (bundle/'parts').glob('*.obj')})
        files.update({str(p.relative_to(package)):p for p in package.rglob('*.so*') if p.is_file()})
        def job(name,argv,identity,environment=None):
            event('JOB_STARTED',name=name)
            r=run_source_job(out/name,argv,identity,environment=environment or {})
            event('JOB_FINISHED',name=name,result=r.status,returncode=r.returncode)
            return r.status=='PROCESS_EXITED' and r.returncode==0
        if not job('physics-job-01',[str(python),'-B',str(HERE/'flower_mujoco.py'),'--bundle',str(bundle),'--output','{job}/output'],files):
            raise RuntimeError('MUJOCO_MATRIX_FAILED')
        physics=out/'physics-job-01/output'
        renderfiles={'renderer':HERE/'flower_mujoco_blender.py','blender':Path('/home/lin/.local/opt/blender-5.1.2-linux-x64/blender'),**files}
        renderfiles.update({p.name:p for p in physics.iterdir() if p.suffix in ('.npz','.json')})
        rendered=job('render-job-01',['/home/lin/.local/bin/blender','-b','-t','4','--python',str(HERE/'flower_mujoco_blender.py'),'--','--bundle',str(bundle),'--physics',str(physics),'--output','{job}/output'],renderfiles)
        if rendered:
            render=out/'render-job-01/output'
            font='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
            vf=f"drawtext=fontfile={font}:text='MuJoCo真实轨迹 · Blender渲染 · 优化前 / 优化后':x=20:y=18:fontsize=23:fontcolor=white:box=1:boxcolor=black@0.45,drawtext=fontfile={font}:text='正常速度':x=20:y=555:fontsize=22:fontcolor=white:enable='lt(t,8)',drawtext=fontfile={font}:text='0.05倍慢放 · 位移未放大':x=20:y=555:fontsize=22:fontcolor=white:enable='gte(t,8)'"
            encodefiles={p.name:p for p in (render/'frames').glob('*.png')};encodefiles['render_manifest']=render/'render.json'
            if job('video-job-01',['/usr/bin/ffmpeg','-nostdin','-n','-framerate','24','-i',str(render/'frames/frame-%04d.png'),'-vf',vf,'-c:v','libx264','-threads','4','-crf','19','-pix_fmt','yuv420p','{job}/flower-mujoco.mp4'],encodefiles):
                html='''<!doctype html><meta charset="utf-8"><title>MuJoCo + Blender 花盘</title><style>body{font:18px/1.7 system-ui;max-width:1100px;margin:auto;padding:30px;background:#f2f0ea}img,video{max-width:100%}a{color:#24547a}</style><h1>同一花盘 · MuJoCo + Blender</h1><p>重建方法验证。64场景：理想轴承／10 Hz、阻尼比0.1柔性支撑，八相位、两个时间步、优化前后。无持续驱动、无额外轴承摩擦；非实测标定。</p><video controls src="video-job-01/flower-mujoco.mp4" poster="render-job-01/output/preview.png"></video><p>真实柔性支撑轨迹，首相位；前8秒正常速度，随后0.05倍慢放。位移没有人为放大。桌面上的支架承载花盘，不是尖底陀螺的自由桌面接触。</p><h2>三层展开（仅展示）</h2><img src="render-job-01/output/three-layer-exploded-display-only.png"><p>上下灰色盖层与中间材料均存在；展开变换不进入物理计算。</p><p><a href="render-job-01/output/flower-mujoco-replay.blend">Blender场景</a> · <a href="physics-job-01/output/summary.json">全部物理统计</a> · <a href="mujoco-results.zip">轨迹与MJCF包</a> · <a href="../table-scene-v1/flower-complete-obj.zip">完整OBJ</a></p>'''
                (out/'index.html').write_text(html)
                with zipfile.ZipFile(out/'mujoco-results.zip','x',zipfile.ZIP_DEFLATED) as z:
                    for p in physics.iterdir():z.write(p,'physics/'+p.name)
                    for p in (HERE/'flower_mujoco.py',HERE/'flexible_bearing.py'):z.write(p,p.name)
                    z.write(bundle/'table-scene.json','table-scene.json')
                    for p in (bundle/'parts').glob('*.obj'):z.write(p,'parts/'+p.name)
                event('FLOWER_RENDER_AND_PHYSICS_PUBLISHED',page=str(out/'index.html'))
        # Preserved failed boundary is inspected only, never promoted to a stage.
        native=ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        boundary=HERE/'noodles-dense-v5-domain-stage-04/output/design.rational.json'
        job('noodles-diagnostic-job-01',['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(HERE/'profile_noodles_final_region.py'),'--source',str(boundary),'--output','{job}/output'],
            {'worker':HERE/'profile_noodles_final_region.py','boundary':boundary,'native':native,'snapshot':ROOT/'python/make_it_stand/verification/spin_problem.py'},
            {'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    finally:
        signal_coordinator(signal.SIGCONT);event('HEART_COORDINATOR_RESUMED')


if __name__=='__main__':main()
