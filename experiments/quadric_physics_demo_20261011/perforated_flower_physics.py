"""Saved three-layer method result, independent moments and ideal-bearing replay.

No contact or manufacturing claim: display arcs are polygonized separately.
"""
from pathlib import Path
import json
import sys

def original_entities(problem):
    layers = [dict(loops=problem['design_region_loops'],
                   z_interval=problem['extrusion_interval'],
                   volume_density=problem['extrusion_volume_density']),
              *problem['fixed_extrusion_layers']]
    return [dict(geometry_representation='exact_polygonal_prism', **x) for x in layers]

def main(folder, out):
    import hashlib
    import subprocess
    import numpy as np
    import pybullet as p
    from PIL import Image, ImageDraw, ImageFont
    from scipy.spatial.transform import Rotation
    here = Path(__file__).resolve().parent
    sys.path[:0] = [str(here), str(here.parents[1]/'python')]
    from bearing_physics import simulate
    from solid_view import extrude, entity_loops
    from make_it_stand.verification.planar_material import integrate_planar_entity
    applications = [json.loads((folder/f'run-{i}/application.json').read_text()) for i in (1,2)]
    case = applications[0]['cases'][0]
    if case['scientific_hash'] != applications[1]['cases'][0]['scientific_hash']:
        raise ValueError('INDEPENDENT_REPEAT_MISMATCH')
    for i in (1,2):
        audit = json.loads((folder/f'check-{i}/scientific.json').read_text())
        if not all(a['moments_match'] and a['balance_constraints_pass'] for a in audit['cases'][0]['attempts']):
            raise ValueError('INDEPENDENT_MOMENT_OR_BALANCE_FAILURE')
    attempt = next(a for a in case['attempts'] if a['seed']==0 and a['status']=='CONVERGED')
    if attempt['solver_result']['accepted_steps'] <= 0:
        raise ValueError('NO_ACCEPTED_UPDATES')
    entities = dict(original=original_entities(case['planar_problem']),
                    optimized=attempt['solver_result']['extruded_material_entities'])
    out.mkdir(parents=True, exist_ok=False)
    def save(name, value):
        (out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    props, rows, files = {}, {}, {}
    for name, parts in entities.items():
        raw = sum((np.array(integrate_planar_entity(e,origin=[0,0],length_scale=3)['raw_mass_moments'])
                   for e in parts), start=np.zeros(10))
        v = raw.copy(); v[0]*=1200*.02**3; v[1:4]*=1200*.02**4; v[4:]*=1200*.02**5
        m=v[0]; c=v[1:4]/m
        Q=np.array([[v[4],v[9],v[8]],[v[9],v[5],v[7]],[v[8],v[7],v[6]]])
        I=np.trace(Q)*np.eye(3)-Q-m*(np.dot(c,c)*np.eye(3)-np.outer(c,c))
        props[name]=dict(mass_kg=m,center_m=c.tolist(),inertia_kg_m2=I.tolist())
        for dt in (.0005,.00025):
            row=simulate(m,c,I,dt=dt,duration=8,hz=35)
            forces=np.array([r['reaction'][:2] for r in row['records'][1:]])
            row['mean_horizontal_bearing_force_N']=float(np.linalg.norm(forces,axis=1).mean())
            row['analytic_horizontal_force_N']=float(m*np.linalg.norm(c[:2])*(2*np.pi*35)**2)
            if row['maximum_relative_energy_drift']>1e-8 or abs(row['mean_horizontal_bearing_force_N']-row['analytic_horizontal_force_N'])>1e-5:
                raise ValueError('BEARING_PHYSICS_CHECK_FAILED')
            save(f'{name}-{dt}.json',row)
            if dt==.00025: rows[name]=row
        files[name]=[]
        for j,e in enumerate(parts):
            mesh=extrude(entity_loops(e),*e['z_interval']); mesh.vertices*=.02
            path=out/f'{name}-{j}-display.obj'; mesh.export(path,digits=17); files[name].append(path)
    save('properties.json',props)
    save('summary.json',dict(classification='RECONSTRUCTED_METHOD_VALIDATION',
        selected_seed=0, selection='PREREGISTERED_NOT_BEST_RUN',
        identical_scientific_repeat=True,
        source_sha256=hashlib.sha256((folder/'run-1/application.json').read_bytes()).hexdigest(),
        scope='IDEAL_BEARING_PHYSICS_NO_CONTACT_NO_MANUFACTURING_APPROVAL',
        display='POLYGONIZED_ARCS_DISPLAY_ONLY',
        outcomes={k:{a:b for a,b in v.items() if a!='records'} for k,v in rows.items()}))
    client=p.connect(p.DIRECT)
    bodies={}
    for name, paths in files.items():
        bodies[name]=[]
        for j,path in enumerate(paths):
            color=([.92,.60,.21,1] if name=='original' else [.05,.66,.59,1]) if j==0 else [.67,.72,.79,1]
            visual=p.createVisualShape(p.GEOM_MESH,fileName=str(path),rgbaColor=color)
            bodies[name].append(p.createMultiBody(baseMass=0,baseVisualShapeIndex=visual,baseCollisionShapeIndex=-1,basePosition=[10,0,0]))
    view=p.computeViewMatrix([.059,-.069,.064],[0,0,0],[0,0,1])
    proj=p.computeProjectionMatrixFOV(39,480/360,.001,2)
    font=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',20)
    small=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',16)
    with (out/'video.log').open('wb') as log:
        proc=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24','-video_size','960x490','-framerate','30','-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','20','-pix_fmt','yuv420p',str(out/'flower-comparison.mp4')],stdin=subprocess.PIPE,stdout=log,stderr=log)
        for frame in range(541):
            slow=frame>=240; t=(frame-240)/300 if slow else frame/30
            canvas=Image.new('RGB',(960,490),'#edf2f8'); draw=ImageDraw.Draw(canvas)
            for j,name in enumerate(rows):
                for ids in bodies.values():
                    for body in ids:p.resetBasePositionAndOrientation(body,[10,0,0],[0,0,0,1])
                state=rows[name]['records'][min(round(t/.00025),32000)]
                quat=Rotation.from_rotvec([0,0,state['angle']]).as_quat()
                for body in bodies[name]:p.resetBasePositionAndOrientation(body,[0,0,0],quat)
                pixels=p.getCameraImage(480,360,view,proj,renderer=p.ER_TINY_RENDERER,lightDirection=[-2,-3,5])[2]
                canvas.paste(Image.fromarray(np.array(pixels,np.uint8).reshape(360,480,4)[:,:,:3]),(480*j,50))
                draw.text((480*j+15,12),'原始三层镂空花盘' if j==0 else 'Faster 优化后的三层花盘',font=font,fill='#18324b')
                draw.text((480*j+15,418),f'轴承水平载荷 {rows[name]["mean_horizontal_bearing_force_N"]:.4g} N',font=small,fill='#18324b')
            draw.text((15,447),f'PyBullet 实际轨迹 / t={t:.3f}s / '+('0.1倍慢放' if slow else '35Hz正常速度，可能有采样错觉'),font=small,fill='#18324b')
            draw.text((15,471),'重建方法案例；理想固定轴承，无驱动、无阻尼；非接触或停止时间预测。',font=small,fill='#18324b')
            if frame==0: canvas.save(out/'preview.png')
            proc.stdin.write(np.asarray(canvas).tobytes())
        proc.stdin.close()
        if proc.wait(timeout=60):raise RuntimeError('VIDEO_ENCODING_FAILED')
    p.disconnect(client)
    print(json.dumps(props),flush=True)

if __name__=='__main__':
    main(Path(sys.argv[1]),Path(sys.argv[2]))
