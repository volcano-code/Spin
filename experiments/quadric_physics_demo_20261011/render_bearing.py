"""Saved planar material -> independent SI moments -> Bullet -> video."""
from pathlib import Path
import sys,json,hashlib,subprocess,io
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path[:0]=[str(HERE),str(ROOT/'python')]
from bearing_physics import simulate
from make_it_stand.verification.planar_material import integrate_planar_entity

def main():
    from PIL import Image,ImageDraw,ImageFont
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MP
    from matplotlib.patches import PathPatch
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    folder=HERE/'multiple-cases-v1';source=folder/'flower-run-1/application.json'
    app=json.loads(source.read_text());s=app['cases'][0]['attempts'][0]['solver_result']
    g=json.loads((folder/'flower-geometry.json').read_text())
    original=dict(geometry_representation='exact_polygonal_prism',loops=[g['outer'],g['bore'][::-1]],z_interval=[-.075,.075],volume_density=1.)
    entities={'original':[original],'optimized':s['extruded_material_entities']}
    props={};results=[]
    protocol=dict(duration_s=8,initial_spin_hz=35,timesteps=[.0005,.00025],
        bearing='ideal revolute at origin / vertical z',gravity=9.81,drive=False,damping=0,
        contact=False,material_scope='INDEPENDENT_MOMENTS_NOT_FULL_TOPOLOGY_OR_MANUFACTURING',
        input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        meters_per_unit=.020,density_kg_m3=1200,video='actual Bullet angle, 0.05x slow playback; 160 seconds for 8 physical seconds')
    for name,parts in entities.items():
        raw=sum((np.array(integrate_planar_entity(e,origin=[0,0],length_scale=3)['raw_mass_moments']) for e in parts),start=np.zeros(10))
        # Independent signed boundary integration, then dimensional SI map.
        v=raw.copy();v[0]*=1200*.02**3;v[1:4]*=1200*.02**4;v[4:]*=1200*.02**5
        m=v[0];c=v[1:4]/m
        Q=np.array([[v[4],v[9],v[8]],[v[9],v[5],v[7]],[v[8],v[7],v[6]]])
        I=np.trace(Q)*np.eye(3)-Q-m*(np.dot(c,c)*np.eye(3)-np.outer(c,c))
        props[name]=dict(mass=m,center=c.tolist(),inertia=I.tolist(),raw_moments=v.tolist())
        for dt in protocol['timesteps']:
            row=simulate(m,c,I,dt=dt,duration=8,hz=35);row['material']=name
            # Static vertical weight remains; compare horizontal bearing load.
            forces=np.array([x['reaction'][:2] for x in row['records'][1:]])
            row['mean_horizontal_bearing_force_N']=float(np.linalg.norm(forces,axis=1).mean())
            row['analytic_horizontal_force_N']=float(m*np.linalg.norm(c[:2])*(2*np.pi*35)**2)
            row['maximum_speed_error']=float(max(abs(x['omega']-2*np.pi*35) for x in row['records']))
            (out/f'{name}-{dt}.json').write_text(json.dumps(row));results.append(row)
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2));(out/'properties.json').write_text(json.dumps(props,indent=2))
    summary=[{k:v for k,v in r.items() if k!='records'} for r in results]
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    for r in results:
        assert r['maximum_relative_energy_drift']<1e-8
        assert abs(r['mean_horizontal_bearing_force_N']-r['analytic_horizontal_force_N'])<1e-5
    def loops(e):
        if 'loops' in e:return e['loops']
        answer=[]
        for cycle in e['boundary_cycles']:
            pts=[]
            for seg in cycle:
                if seg['type']=='LINE_SEGMENT':pts.append(seg['start'])
                else:
                    t=seg['start_angle']+np.linspace(0,seg['sweep'],max(4,int(abs(seg['sweep'])*100)),endpoint=False)
                    pts.extend(np.array(seg['center'])+seg['radius']*np.column_stack([np.cos(t),np.sin(t)]))
            answer.append(pts)
        return answer
    sprites={}
    for name,parts in entities.items():
        fig=plt.figure(figsize=(4,4),dpi=120);ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(-1.5,1.5),ylim=(-1.5,1.5));ax.axis('off');ax.set_aspect('equal')
        for j,e in enumerate(parts):
            verts=[];codes=[]
            for loop in loops(e):
                q=np.asarray(loop);verts.extend([*q,q[0]]);codes.extend([MP.MOVETO]+[MP.LINETO]*(len(q)-1)+[MP.CLOSEPOLY])
            ax.add_patch(PathPatch(MP(verts,codes),facecolor=('#e2a953' if name=='original' else '#26a99d') if j==0 else '#9baabe',edgecolor='#24364c',lw=.8))
        buffer=io.BytesIO();fig.savefig(buffer,format='png',transparent=True);plt.close(fig);buffer.seek(0);sprites[name]=Image.open(buffer).convert('RGBA')
    selected=[r for r in results if r['dt']==.00025]
    font=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',21)
    small=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',17)
    # Slow motion samples actual trajectory densely enough to avoid 35Hz/30fps aliasing.
    with (out/'ffmpeg.log').open('wb') as log:
        proc=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24','-video_size','960x610','-framerate','30','-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','21','-pix_fmt','yuv420p',str(out/'flower-bearing.mp4')],stdin=subprocess.PIPE,stdout=log,stderr=log)
        for frame in range(4801):
            t=min(frame/600,8);canvas=Image.new('RGB',(960,610),'#edf2f8');draw=ImageDraw.Draw(canvas)
            for j,row in enumerate(selected):
                state=row['records'][min(round(t/.00025),len(row['records'])-1)]
                sprite=sprites[row['material']].rotate(np.rad2deg(state['angle'])%360,resample=Image.Resampling.BICUBIC)
                canvas.paste(sprite,(480*j,65),sprite)
                draw.text((480*j+18,12),'原始偏心五瓣' if j==0 else 'Faster 二维优化后',font=font,fill='#17334c')
                draw.text((480*j+18,45),f'质量 {row["mass"]*1000:.3f} g  |  转速 {state["omega"]/(2*np.pi):.2f} Hz',font=small,fill='#17334c')
                x,y=480*j+240,305;draw.ellipse((x-6,y-6,x+6,y+6),fill='#24364c')
                draw.text((480*j+18,520),f'轴承水平载荷 {row["mean_horizontal_bearing_force_N"]:.6g} N',font=small,fill='#17334c')
            draw.text((18,554),f'物理时间 {t:.3f} / 8 s · 0.05倍慢放 · PyBullet 实际关节轨迹',font=small,fill='#17334c')
            draw.text((18,583),'理想轴承 / 无驱动、无阻尼；两者不停止。效果看轴承载荷，不看固定轴的摆动。',font=small,fill='#17334c')
            if frame==0:canvas.save(out/'preview.png')
            proc.stdin.write(np.asarray(canvas).tobytes())
        proc.stdin.close();assert proc.wait(timeout=60)==0
    assert hashlib.sha256(source.read_bytes()).hexdigest()==protocol['input_sha256']
    print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
