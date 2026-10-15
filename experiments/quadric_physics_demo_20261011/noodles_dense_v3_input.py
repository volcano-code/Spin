"""Separate preregistered v3: dense rounded tendrils, not optimized material."""
from pathlib import Path
import hashlib
import json
import numpy as np
import trimesh
from noodles_input import tube_mesh,mesh_audit,write_obj

REFERENCE=Path('/tmp/codex-clipboard-badb96ea-8a95-436b-9dcc-7aae4b444789.png')
SPEC=dict(schema='dense-noodles-reconstruction-v3',classification='RECONSTRUCTED_METHOD_VALIDATION',
    author_original=False,source_view='USER_REFERENCE_LEFT_UNOPTIMIZED_ONLY',
    paper_locator='Spin-It Faster Section 6.2; reconstructed visible arrangement only',
    ellipsoid_center_model=[0.,0.,.90],ellipsoid_radii_model=[.98,.82,.70],inner_scale=.94,
    ellipsoid_subdivisions=3,tube_radius_model=.030,tube_samples=33,tube_sections=12,tube_cap_rings=4,
    tube_z_rows=11,tube_z_span_model=.57/.65*.70,tube_y_lanes=9,tube_y_span_model=.66/.80*.82,
    layout_ellipsoid_norm_limit=.88,tube_end_extension_model=.035,tube_wave_z_model=.035,
    tube_wave_y_model=.014,tube_wave_phase_increment=.37,tube_surface_arch_y_model=.60,
    hub_radius_model=.056,hub_interval_model=[.21,1.59],shaft_radius_model=.045,
    grip_radius_model=.105,grip_center_z_model=1.78,tip_bottom_model=-.16,
    model_unit_m=.025,density_kg_m3=1200.,target_origin_model=[0.,0.,-.16],target_axis=[0.,0.,1.],
    objective='FASTER_INERTIA_RATIO_PLUS_CENTER_OF_MASS_HEIGHT',center_of_mass_weight=1.,
    weight_coordinate_contract='FASTER_NORMALIZED_REQUEST',backend='CGAL',normal_subdivision=2,seeds=[0,1],
    maximum_attempts_per_seed=8,maximum_evaluations_per_attempt=60,raw_residual_tolerance=1e-12,alpha=.5,t_max=.5,
    parameters_source='RECONSTRUCTED_HIDDEN_GEOMETRY_FIXED_BEFORE_OPTIMIZATION',
    hidden_geometry_assumptions=['rear tendril layout','internal connections','curves and rounded tube caps',
        'ellipsoid candidate dimensions','axle and grip dimensions'],
    visual_change_purpose='Rounder envelope and visible horizontal tendrils; concealed thinner hub; no cut preference')


def tube_layout():
    rows=[];radii=np.asarray(SPEC['ellipsoid_radii_model'])
    for z in np.linspace(-SPEC['tube_z_span_model'],SPEC['tube_z_span_model'],SPEC['tube_z_rows']):
        for y in np.linspace(-SPEC['tube_y_span_model'],SPEC['tube_y_span_model'],SPEC['tube_y_lanes']):
            level=(y/radii[1])**2+(z/radii[2])**2
            if level<SPEC['layout_ellipsoid_norm_limit']:
                rows.append(dict(y_model=float(y),z_model=float(SPEC['ellipsoid_center_model'][2]+z),
                    half_length_model=float(radii[0]*np.sqrt(1-level)+SPEC['tube_end_extension_model'])))
    return rows


def dense_centerline(index):
    row=tube_layout()[index];u=np.linspace(0.,1.,SPEC['tube_samples'])
    phase=SPEC['tube_wave_phase_increment']*index
    x=(2*u-1)*row['half_length_model']
    y=(row['y_model']*(2*u-1)**2+
       SPEC['tube_surface_arch_y_model']*row['y_model']/SPEC['tube_y_span_model']*np.sin(2*np.pi*u)**2+
       SPEC['tube_wave_y_model']*np.sin(4*np.pi*u)*np.sin(np.pi*u))
    z=row['z_model']+SPEC['tube_wave_z_model']*np.sin(3*np.pi*u+phase)*np.sin(np.pi*u)
    return np.column_stack((x,y,z))


def identity():
    return hashlib.sha256(json.dumps(dict(spec=SPEC,tube_layout=tube_layout()),sort_keys=True,
        separators=(',',':'),allow_nan=False).encode()).hexdigest()


def rounded_tube(points,radius,sections):
    """Reuse parallel-transport body, replace disk caps with round hemispheres."""
    base=tube_mesh(points,radius,sections)
    rings=base.vertices[:-2].reshape(len(points),sections,3)
    tangent=np.gradient(points,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    theta=np.linspace(0,np.pi/2,SPEC['tube_cap_rings']+2)[1:-1]
    first=[points[0]-radius*np.sin(t)*tangent[0]+np.cos(t)*(rings[0]-points[0]) for t in theta[::-1]]
    last=[points[-1]+radius*np.sin(t)*tangent[-1]+np.cos(t)*(rings[-1]-points[-1]) for t in theta]
    all_rings=np.asarray(first+list(rings)+last)
    vertices=np.vstack((all_rings.reshape(-1,3),points[0]-radius*tangent[0],points[-1]+radius*tangent[-1]))
    faces=[]
    for i in range(len(all_rings)-1):
        for j in range(sections):
            a=i*sections+j;b=i*sections+(j+1)%sections;c=(i+1)*sections+j;d=(i+1)*sections+(j+1)%sections
            faces.extend(((a,b,d),(a,d,c)))
    start,end=len(vertices)-2,len(vertices)-1
    for j in range(sections):
        k=(j+1)%sections;last_ring=(len(all_rings)-1)*sections
        faces.extend(((start,k,j),(end,last_ring+j,last_ring+k)))
    mesh=trimesh.Trimesh(vertices,np.asarray(faces),process=False)
    if mesh.volume<0:mesh.invert()
    return mesh


def components():
    sphere=trimesh.creation.icosphere(subdivisions=SPEC['ellipsoid_subdivisions'])
    center=np.asarray(SPEC['ellipsoid_center_model'])
    sphere.vertices=sphere.vertices*np.asarray(SPEC['ellipsoid_radii_model'])+center
    inner=sphere.copy();inner.vertices=center+SPEC['inner_scale']*(inner.vertices-center)
    hub=trimesh.creation.cylinder(radius=SPEC['hub_radius_model'],height=np.diff(SPEC['hub_interval_model'])[0],sections=24)
    hub.apply_translation([0,0,np.mean(SPEC['hub_interval_model'])])
    profile=[(0.,-.16)]+[(.025*np.cos(a),-.135+.025*np.sin(a)) for a in np.linspace(-np.pi/2,0,7)[1:]]
    r,z=SPEC['grip_radius_model'],SPEC['grip_center_z_model'];start=-np.arccos(SPEC['shaft_radius_model']/r)
    profile +=[(SPEC['shaft_radius_model'],-.03),(SPEC['shaft_radius_model'],z+r*np.sin(start))]
    profile +=[(r*np.cos(a),z+r*np.sin(a)) for a in np.linspace(start,np.pi/2,17)[1:-1]]+[(0.,z+r)]
    pin=trimesh.creation.revolve(np.asarray(profile),sections=32)
    if pin.volume<0:pin.invert()
    parts={'outer-ellipsoid':sphere,'inner-ellipsoid-candidate':inner,'hub-fixed':hub,'pin-fixed':pin}
    for i in range(len(tube_layout())):
        parts[f'tube-{i:03d}-fixed']=rounded_tube(dense_centerline(i),SPEC['tube_radius_model'],SPEC['tube_sections'])
    return parts


def preview(directory,parts):
    import pybullet as p
    from PIL import Image,ImageDraw,ImageFont
    canvas=Image.new('RGB',(1600,1100),'#f4f5f7');draw=ImageDraw.Draw(canvas)
    font=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',23)
    small=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',18)
    draw.text((20,12),'面条 v3｜未优化初始材料：圆润管束＋椭球候选＋握柄与底尖',font=font,fill='#243249')
    client=p.connect(p.DIRECT)
    try:
        for name in parts:
            if name=='inner-ellipsoid-candidate':continue
            color=[.75,.27,.49,1.] if name.startswith('tube-') else [.86,.47,.65,1.]
            if name=='pin-fixed':color=[.39,.42,.46,1.]
            shape=p.createVisualShape(p.GEOM_MESH,fileName=str(directory/(name+'-model.obj')),rgbaColor=color,physicsClientId=client)
            p.createMultiBody(0,baseVisualShapeIndex=shape,baseCollisionShapeIndex=-1,physicsClientId=client)
        cameras=[([0,-5,.92],[0,0,.88]),([5,0,.92],[0,0,.88]),([3.6,-5,3.0],[0,0,.88])]
        labels=['正面：实际初始材料包含椭球','侧面：端部圆帽与前后弯曲','透视：内部细连接隐藏在初始实体中']
        projection=p.computeProjectionMatrixFOV(32,800/440,.01,20)
        for i,(eye,target) in enumerate(cameras):
            pixels=p.getCameraImage(800,440,p.computeViewMatrix(eye,target,[0,0,1]),projection,
                renderer=p.ER_TINY_RENDERER,lightDirection=[-2,-3,5],physicsClientId=client)[2]
            panel=Image.fromarray(np.asarray(pixels,np.uint8).reshape(440,800,4)[:,:,:3])
            x=(i%2)*800;y=90+(i//2)*480
            canvas.paste(panel,(x,y));draw.text((x+20,y-30),labels[i],font=font,fill='#243249')
    finally:p.disconnect(client)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.4),dpi=100)
    for name,mesh in parts.items():
        section=mesh.section(plane_origin=[0,0,0],plane_normal=[0,1,0])
        if section is None:continue
        for line in section.discrete:
            if name=='outer-ellipsoid':ax.fill(line[:,0],line[:,2],color='#dea0ba',alpha=.5)
            elif name=='inner-ellipsoid-candidate':ax.plot(line[:,0],line[:,2],'--',color='#b05282',linewidth=1.3)
            else:ax.plot(line[:,0],line[:,2],color='#713854' if name.startswith('tube-') else '#41454b',linewidth=.7)
    ax.set(xlim=(-1.15,1.15),ylim=(-.2,1.95),xlabel='x [model]',ylabel='z [model]')
    ax.set_aspect('equal');ax.set_title('Actual y=0 source sections; dashed: candidate domain')
    fig.tight_layout();fig.canvas.draw();pixels=np.asarray(fig.canvas.buffer_rgba())[:,:,:3];plt.close(fig)
    canvas.paste(Image.fromarray(pixels).resize((800,440)),(800,570))
    draw.text((820,540),'实际网格剖面：未预制任何优化空腔',font=font,fill='#243249')
    draw.text((20,1035),'79 根圆帽细管；隐藏形状为预登记重建假设，非作者原始资产。',font=small,fill='#243249')
    draw.text((20,1065),'各部件已做初检；精确整体并集、分区、可达性、优化及物理仿真尚未通过。',font=small,fill='#243249')
    canvas.save(directory/'input-preview.png')


def prepare(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    try:
        parts=components();entries={}
        for name,mesh in parts.items():
            audit=mesh_audit(mesh)
            if not audit['closed'] or not audit['winding_consistent'] or audit['components']!=1 or mesh.volume<=0 or audit['minimum_triangle_area_model2']<=0:
                raise ValueError('INVALID_SOURCE_PART: '+name)
            write_obj(directory/(name+'-model.obj'),mesh)
            si=mesh.copy();si.vertices*=SPEC['model_unit_m'];write_obj(directory/(name+'-si.obj'),si)
            entries[name]=dict(audit=audit,model_obj=name+'-model.obj',si_obj=name+'-si.obj',
                sha256={s:hashlib.sha256((directory/(name+'-'+s+'.obj')).read_bytes()).hexdigest() for s in ('model','si')})
        manifest=dict(schema='dense-noodles-source-parts-v3',spec=SPEC,spec_sha256=identity(),tube_layout=tube_layout(),
            reference_image_sha256=hashlib.sha256(REFERENCE.read_bytes()).hexdigest() if REFERENCE.exists() else None,
            parts=entries,construction=dict(fixed_parts=['hub-fixed','pin-fixed']+[n for n in parts if n.startswith('tube-')],
                outer_operation='UNION(outer-ellipsoid, ALL fixed_parts)',
                design_operation='DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',
                fixed_material_operation='DIFFERENCE(outer, design)'),
            connections='TENDRILS_PENETRATE_INTERNAL_HUB; EXACT_WHOLE_UNION_PENDING',
            preview_role='ALL_INITIAL_MATERIAL_COMPONENTS_VISIBLE_NOT_OPTIMIZED_BOUNDARY',
            status='SOURCE_PARTS_PREPARED_PENDING_EXACT_ADMISSION',optimization_invocations=0,simulation_invocations=0,
            creator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            tube_mesh_helper_sha256=hashlib.sha256(Path(__file__).with_name('noodles_input.py').read_bytes()).hexdigest())
        (directory/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
        preview(directory,parts)
        return manifest
    except Exception as exc:
        (directory/'failure.json').write_text(json.dumps(dict(status='SOURCE_PREPARATION_FAILED',
            type=type(exc).__name__,message=str(exc),optimization_invocations=0,simulation_invocations=0),indent=2)+'\n')
        raise


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',type=Path,required=True)
    prepare(parser.parse_args().prepare)
