"""Real 3D flower trajectory replay and a newly built solid handled Heart."""
from pathlib import Path
import sys,json,subprocess,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path[:0]=[str(HERE),str(ROOT/'python'),str(ROOT/'experiments/faster_heart_contact_20261011'),str(ROOT/'experiments/heart_stem_contact')]
from solid_view import extrude,entity_loops

def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False))

def main():
    import trimesh,pybullet as p
    from scipy.spatial.transform import Rotation
    from PIL import Image,ImageDraw,ImageFont
    from rounded_input import rounded_components,SPEC
    from stem_geometry import boolean
    from demo import to_physics
    from bullet_contact_common_frame import create_body
    from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    source=HERE/'multiple-cases-v1';app=json.loads((source/'flower-run-1/application.json').read_text())
    s=app['cases'][0]['attempts'][0]['solver_result'];g=json.loads((source/'flower-geometry.json').read_text())
    entities={'original':[dict(loops=[g['outer'],g['bore'][::-1]],z_interval=[-.075,.075])],
              'optimized':s['extruded_material_entities']}
    files={};audits={}
    for name,parts in entities.items():
        files[name]=[];audits[name]=[]
        for j,e in enumerate(parts):
            mesh=extrude(entity_loops(e),*e['z_interval']);mesh.vertices*=.02
            path=out/f'flower-{name}-{j}.obj';mesh.export(path,digits=17);files[name].append(path)
            audits[name].append(dict(vertices=len(mesh.vertices),triangles=len(mesh.faces),closed=mesh.is_watertight,volume_m3=mesh.volume,
                components=len(mesh.split()),role='DISPLAY_DERIVATIVE_ARCS_POLYGONIZED_NOT_AUTHORITY_OR_COLLISION'))
    save(out/'flower-display-audit.json',audits)
    client=p.connect(p.DIRECT)
    font=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',21)
    small=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',16)
    def meshbody(file,color):
        shape=p.createVisualShape(p.GEOM_MESH,fileName=str(file),rgbaColor=color,physicsClientId=client)
        return p.createMultiBody(baseMass=0,baseVisualShapeIndex=shape,baseCollisionShapeIndex=-1,physicsClientId=client)
    flower={}
    for name,paths in files.items():
        flower[name]=[meshbody(f,([.9,.56,.16,1] if name=='original' else [.03,.63,.58,1]) if j==0 else [.55,.61,.7,1]) for j,f in enumerate(paths)]
    # Fixed bearing support is a visual representation of the ideal joint,
    # not extra mass or a second collision model.
    visual=p.createVisualShape(p.GEOM_CYLINDER,radius=.0018,length=.016,rgbaColor=[.2,.24,.3,1],physicsClientId=client)
    shaft=p.createMultiBody(baseMass=0,baseVisualShapeIndex=visual,basePosition=[0,0,-.008],physicsClientId=client)
    basevisual=p.createVisualShape(p.GEOM_CYLINDER,radius=.009,length=.004,rgbaColor=[.23,.28,.36,1],physicsClientId=client)
    stand=p.createMultiBody(baseMass=0,baseVisualShapeIndex=basevisual,basePosition=[0,0,-.016],physicsClientId=client)
    rows={name:json.loads((source/f'flower-bearing-job-v2/output/{name}-0.00025.json').read_text()) for name in flower}
    view=p.computeViewMatrix([.066,-.075,.072],[0,0,-.001],[0,0,1]);proj=p.computeProjectionMatrixFOV(39,480/360,.001,2)
    def image():
        a=p.getCameraImage(480,360,viewMatrix=view,projectionMatrix=proj,renderer=p.ER_TINY_RENDERER,
            lightDirection=[-2,-3,5],shadow=1,physicsClientId=client)[2]
        return Image.fromarray(np.asarray(a,np.uint8).reshape(360,480,4)[:,:,:3])
    def encode(path,frames):
        with path.with_suffix('.log').open('wb') as log:
            proc=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24','-video_size','960x490','-framerate','30','-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','20','-pix_fmt','yuv420p',str(path)],stdin=subprocess.PIPE,stdout=log,stderr=log)
            for im in frames:proc.stdin.write(np.asarray(im).tobytes())
            proc.stdin.close()
            if proc.wait(timeout=60):raise RuntimeError('VIDEO_ENCODING_FAILED')
    def flower_frames():
        for n in range(541):
            slow=n>=240;t=n/30 if not slow else (n/30-8)*.1
            canvas=Image.new('RGB',(960,490),'#edf2f8');d=ImageDraw.Draw(canvas)
            for j,name in enumerate(flower):
                for ids in flower.values():
                    for id in ids:p.resetBasePositionAndOrientation(id,[10,0,0],[0,0,0,1],physicsClientId=client)
                r=rows[name]['records'];k=min(round(t/.00025),len(r)-1);angle=r[k]['angle']
                for id in flower[name]:p.resetBasePositionAndOrientation(id,[0,0,0],Rotation.from_rotvec([0,0,angle]).as_quat(),physicsClientId=client)
                canvas.paste(image(),(480*j,50))
                d.text((480*j+16,12),'原始五瓣 · 三维实体' if j==0 else 'Faster 优化后 · 三维实体',font=font,fill='#18324b')
                d.text((480*j+16,416),f'水平轴承载荷 {rows[name]["mean_horizontal_bearing_force_N"]:.4g} N',font=small,fill='#18324b')
            d.text((16,447),f'实际 PyBullet 轨迹｜物理 t={t:.3f}s｜'+('0.1倍慢放片段' if slow else '正常速度（35Hz可能出现采样错觉）'),font=small,fill='#18324b')
            d.text((16,471),'厚度 3 mm；固定支架表示理想轴承。无电机、无阻尼，不模拟实际停止时间。',font=small,fill='#18324b')
            if n==0:canvas.save(out/'flower-3d.png')
            yield canvas
    encode(out/'flower-3d.mp4',flower_frames());p.disconnect(client)
    print('Flower 3D complete; constructing connected handled Heart',flush=True)
    body,_,pin=rounded_components()
    body.export(out/'heart-body-model.obj',digits=17);pin.export(out/'heart-grip-model.obj',digits=17)
    heart=boolean(body,pin,'union',out)
    if not heart.is_watertight or not heart.is_winding_consistent or len(heart.split())!=1:raise ValueError('HANDLED_HEART_NOT_SINGLE_CLOSED_BODY')
    raw=[float(x) for x in integrate_exact_rational_raw_moments(heart.vertices.tolist(),heart.faces.tolist(),arithmetic_backend='GMPY2_2_2_1_EXACT_V1')]
    props=to_physics(raw);heart.vertices*=.025;heart.export(out/'heart-handled-si.obj',digits=17)
    hull=heart.convex_hull;hull.export(out/'heart-collision-hull-si.obj',digits=17)
    save(out/'heart-input.json',dict(classification='RECONSTRUCTED_SOLID_HANDLED_BASELINE_NOT_OPTIMIZED',construction='actual CGAL union of body and rounded shaft/grip',generator_spec=SPEC,
        mass_kg=props[0],center_m=props[1].tolist(),inertia_kg_m2=props[2].tolist(),closed=True,components=1,
        design_domain=None,manufacturing_approved=False,collision='convex hull; plane support preserved, contact manifold approximation'))
    # New baseline only: no old cavity or old optimized inertia is used.
    results=[]
    for dt in [.0005,.00025]:
        client=p.connect(p.DIRECT);p.setGravity(0,0,-9.81,physicsClientId=client);p.setTimeStep(dt,physicsClientId=client)
        plane=p.createCollisionShape(p.GEOM_PLANE,physicsClientId=client);p.createMultiBody(0,plane,physicsClientId=client)
        R=Rotation.from_rotvec([np.deg2rad(2),0,0]).as_matrix();origin=np.array([0,0,.0002-(heart.vertices@R.T)[:,2].min()])
        id=create_body(p,client,out/'heart-handled-si.obj',[out/'heart-collision-hull-si.obj'],props,R,origin,35)
        rec=[]
        for k in range(round(8/dt)+1):
            if k:p.stepSimulation(physicsClientId=client)
            if k%max(1,round(1/(240*dt)))==0:
                pos,q=p.getBasePositionAndOrientation(id,physicsClientId=client);v,w=p.getBaseVelocity(id,physicsClientId=client)
                axis=(Rotation.from_quat(q).as_matrix()@props[4].T)[:,2]
                rec.append(dict(t=k*dt,position=list(pos),quaternion=list(q),tilt=float(np.rad2deg(np.arccos(np.clip(axis[2],-1,1)))),omega=list(w)))
        result=dict(dt=dt,records=rec,scope='ORIGINAL_SOLID_HANDLED_HEART_CONTACT_DIAGNOSTIC_NOT_OPTIMIZATION_OR_VALIDATED_STABILITY')
        save(out/f'heart-contact-{dt}.json',result);results.append(result);p.disconnect(client)
    client=p.connect(p.DIRECT)
    id=create_body(p,client,out/'heart-handled-si.obj',[out/'heart-collision-hull-si.obj'],props,np.eye(3),[0,0,.004],0,color=(.74,.16,.38,1))
    planevisual=p.createVisualShape(p.GEOM_BOX,halfExtents=[.2,.2,.001],rgbaColor=[.8,.85,.9,1],physicsClientId=client)
    p.createMultiBody(0,baseVisualShapeIndex=planevisual,basePosition=[0,0,-.001],physicsClientId=client)
    view=p.computeViewMatrix([.105,-.14,.083],[0,0,.025],[0,0,1]);proj=p.computeProjectionMatrixFOV(37,480/360,.001,2)
    def heart_frames():
        for n in range(241):
            t=n/30;canvas=Image.new('RGB',(960,490),'#edf2f8');d=ImageDraw.Draw(canvas)
            for j,row in enumerate(results):
                state=min(row['records'],key=lambda x:abs(x['t']-t))
                p.resetBasePositionAndOrientation(id,state['position'],state['quaternion'],physicsClientId=client)
                canvas.paste(image(),(480*j,50));d.text((480*j+16,12),f'带柄实心基线 · dt={row["dt"]}',font=font,fill='#18324b')
                d.text((480*j+16,415),f'倾角 {state["tilt"]:.2f}°',font=small,fill='#18324b')
            d.text((16,448),f'实际接触运动 t={t:.2f}s｜新带柄完整实体；不是优化后结果',font=small,fill='#18324b')
            d.text((16,472),'35 Hz / 初倾2° / 无驱动 / 单凸包碰撞近似。两时间步并列，不承诺旋转稳定。',font=small,fill='#18324b')
            if n==0:canvas.save(out/'heart-handled.png')
            yield canvas
    encode(out/'heart-handled-contact.mp4',heart_frames());p.disconnect(client)
    save(out/'manifest.json',{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in out.iterdir() if x.is_file()})
    print('3D videos and new solid handled Heart completed',flush=True)

if __name__=='__main__':main()
