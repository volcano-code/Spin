"""Torque-free Bullet comparison of the saved non-principal Cuboid method case.

No string, gravity, pivot, drive, or collision. This tests rotational axis
alignment, not a complete yo-yo toy or real stopping time.
"""
from pathlib import Path
import sys,json,hashlib,subprocess
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path[:0]=[str(HERE),str(ROOT/'experiments/heart_stem_contact')]
from demo import to_physics
from bullet_contact_common_frame import create_body

def main():
    import pybullet as p,trimesh
    from scipy.spatial.transform import Rotation
    from PIL import Image,ImageDraw,ImageFont
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    data=HERE/'multiple-cases-v1';source=data/'cuboid-material.json';m=json.loads(source.read_text())
    props={name:to_physics(m['moments'][key]) for name,key in [('original','outer'),('optimized','material')]}
    mesh=trimesh.load(data/'cuboid-outer.obj',force='mesh',process=False);mesh.vertices*=.025
    assert mesh.is_watertight and mesh.is_winding_consistent
    mesh.export(out/'outer-si.obj',digits=17)
    protocol=dict(duration_s=8,initial_spin_hz=35,dt_seconds=[.0005,.00025],gravity=0,contacts=False,
        fixed_pivot=False,drive=False,damping=0,interpretation='TORQUE_FREE_AXIS_ALIGNMENT_NOT_STRING_YOYO',
        video_dt=.00025,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2));all_results=[]
    for dt in protocol['dt_seconds']:
        for name in props:
            client=p.connect(p.DIRECT);p.setGravity(0,0,0,physicsClientId=client)
            p.setPhysicsEngineParameter(fixedTimeStep=dt,deterministicOverlappingPairs=1,physicsClientId=client)
            pp=props[name];body=create_body(p,client,out/'outer-si.obj',[out/'outer-si.obj'],pp,np.eye(3),-pp[1],35)
            p.setCollisionFilterGroupMask(body,-1,0,0,physicsClientId=client)
            p.resetBaseVelocity(body,linearVelocity=[0,0,0],angularVelocity=[0,0,2*np.pi*35],physicsClientId=client)
            info=p.getDynamicsInfo(body,-1,physicsClientId=client)
            np.testing.assert_allclose(pp[4]@np.diag(info[2])@pp[4].T,pp[2],atol=1e-18,rtol=1e-12)
            records=[];e0=None;l0=None;edrift=0.;ldrift=0.;maximum=0.
            for step in range(round(8/dt)+1):
                if step:p.stepSimulation(physicsClientId=client)
                pos,q=p.getBasePositionAndOrientation(body,physicsClientId=client)
                v,w=p.getBaseVelocity(body,physicsClientId=client);R=Rotation.from_quat(q).as_matrix()
                local=R.T@w;energy=.5*np.dot(pp[3]*local,local);L=R@(pp[3]*local)
                if e0 is None:e0=energy;l0=L.copy()
                edrift=max(edrift,abs(energy/e0-1));ldrift=max(ldrift,np.linalg.norm(L-l0)/np.linalg.norm(l0))
                axis=(R@pp[4].T)[:,2];angle=float(np.rad2deg(np.arccos(np.clip(axis[2],-1,1))));maximum=max(maximum,angle)
                if step%max(1,round(1/(240*dt)))==0 or step==round(8/dt):
                    records.append(dict(t=step*dt,position=list(pos),quaternion=list(q),axis_deviation_degrees=angle,
                        angular_velocity=list(w),energy_J=float(energy),angular_momentum=L.tolist()))
            result=dict(material=name,dt=dt,records=records,maximum_axis_deviation_degrees=maximum,
                relative_energy_drift=float(edrift),relative_world_momentum_drift=float(ldrift),world_constraints=p.getNumConstraints(physicsClientId=client))
            (out/f'{name}-{dt}.json').write_text(json.dumps(result,indent=2));all_results.append(result);p.disconnect(client)
    (out/'summary.json').write_text(json.dumps([{k:v for k,v in x.items() if k!='records'} for x in all_results],indent=2))
    client=p.connect(p.DIRECT);bodies=[]
    for name,color in [('original',(.9,.45,.2,1)),('optimized',(.02,.6,.6,1))]:
        bodies.append(create_body(p,client,out/'outer-si.obj',[out/'outer-si.obj'],props[name],np.eye(3),-props[name][1],35,color))
    chosen=[next(x for x in all_results if x['material']==n and x['dt']==.00025) for n in props]
    fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc';font=ImageFont.truetype(fp,20);small=ImageFont.truetype(fp,16)
    with (out/'ffmpeg.log').open('wb') as log:
        proc=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24','-video_size','960x460','-framerate','30','-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','20','-pix_fmt','yuv420p',str(out/'yoyo-free-rotation.mp4')],stdin=subprocess.PIPE,stdout=log,stderr=log)
        for frame in range(241):
            t=min(frame/30,8);canvas=Image.new('RGB',(960,460),'#edf2f8');draw=ImageDraw.Draw(canvas)
            for j,row in enumerate(chosen):
                state=min(row['records'],key=lambda x:abs(x['t']-t))
                for body in bodies:p.resetBasePositionAndOrientation(body,[10,10,10],[0,0,0,1],physicsClientId=client)
                p.resetBasePositionAndOrientation(bodies[j],state['position'],state['quaternion'],physicsClientId=client)
                view=p.computeViewMatrix([.08,-.1,.08],[0,0,0],[0,0,1]);projection=p.computeProjectionMatrixFOV(40,480/340,.001,3)
                img=p.getCameraImage(480,340,viewMatrix=view,projectionMatrix=projection,renderer=p.ER_TINY_RENDERER,physicsClientId=client)
                canvas.paste(Image.fromarray(np.asarray(img[2],np.uint8).reshape(340,480,4)[:,:,:3]),(480*j,50))
                draw.text((480*j+15,10),'原始非主轴 Cuboid' if j==0 else '二次曲面优化后',font=font,fill='#142944')
                draw.text((480*j+15,399),f'目标轴偏离初始方向 {state["axis_deviation_degrees"]:.3f}°',font=small,fill='#142944')
            draw.text((15,431),f't={t:.2f}s  |  PyBullet 无外力矩旋转；没有绳索、重力或固定轴；dt=0.00025s',font=small,fill='#425975')
            if frame in (0,240):canvas.save(out/f'frame-{frame//30}s.png')
            proc.stdin.write(np.asarray(canvas).tobytes())
        proc.stdin.close();assert proc.wait(timeout=45)==0
    p.disconnect(client)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,4))
    for row in all_results:
        ax.plot([x['t'] for x in row['records']],[x['axis_deviation_degrees'] for x in row['records']],
            '--' if row['dt']==.0005 else '-',label=f'{row["material"]}, dt={row["dt"]}')
    ax.set(xlabel='Time / s',ylabel='Target-axis deviation / degrees',title='Actual torque-free trajectories, both time steps');ax.legend();ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'axis-curves.png',dpi=150);plt.close(fig)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==protocol['source_sha256']
    print(json.dumps([{k:v for k,v in row.items() if k!='records'} for row in all_results]),flush=True)

if __name__=='__main__':main()
