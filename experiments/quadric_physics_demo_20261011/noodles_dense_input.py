"""Dense horizontal curved tendrils reconstructed from the unoptimized view.

The visible reference guides arrangement only. Hidden curves, radii and the
editable ellipsoid are preregistered assumptions, never recovered author CAD.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import trimesh
from noodles_input import tube_mesh, mesh_audit, write_obj

REFERENCE=Path('/tmp/codex-clipboard-badb96ea-8a95-436b-9dcc-7aae4b444789.png')
SPEC=dict(schema='dense-noodles-reconstruction-v2',classification='RECONSTRUCTED_METHOD_VALIDATION',
    author_original=False,source_view='USER_REFERENCE_LEFT_UNOPTIMIZED_ONLY',paper_locator='Spin-It Faster Section 6.2 and Figure 9 row 6',
    ellipsoid_center_model=[0.,0.,1.03],ellipsoid_radii_model=[1.02,.80,.65],inner_scale=.94,
    ellipsoid_subdivisions=3,tube_radius_model=.030,tube_samples=33,tube_sections=12,
    tube_z_rows=11,tube_z_span_model=.57,tube_y_lanes=9,tube_y_span_model=.66,
    layout_ellipsoid_norm_limit=.88,tube_end_extension_model=.065,tube_wave_z_model=.035,
    tube_wave_y_model=.014,tube_wave_phase_increment=.37,hub_radius_model=.09,hub_interval_model=[.32,1.76],
    shaft_radius_model=.052,grip_radius_model=.10,grip_center_z_model=1.94,
    model_unit_m=.025,density_kg_m3=1200.,target_origin_model=[0.,0.,-.16],target_axis=[0.,0.,1.],
    objective='FASTER_INERTIA_RATIO_PLUS_CENTER_OF_MASS_HEIGHT',center_of_mass_weight=1.,
    weight_coordinate_contract='FASTER_NORMALIZED_REQUEST',backend='CGAL',normal_subdivision=2,seeds=[0,1],
    maximum_attempts_per_seed=8,maximum_evaluations_per_attempt=60,raw_residual_tolerance=1e-12,alpha=.5,t_max=.5,
    parameters_source='RECONSTRUCTED_HIDDEN_GEOMETRY_FIXED_BEFORE_OPTIMIZATION',
    hidden_geometry_assumptions=['rear tendril placement','tube centerlines and diameters','central hub and source shell separation',
                                'ellipsoid dimensions','axle and grip dimensions'])


def tube_layout():
    center=np.asarray(SPEC['ellipsoid_center_model']);radii=np.asarray(SPEC['ellipsoid_radii_model'])
    rows=[]
    for z in np.linspace(-SPEC['tube_z_span_model'],SPEC['tube_z_span_model'],SPEC['tube_z_rows']):
        for y in np.linspace(-SPEC['tube_y_span_model'],SPEC['tube_y_span_model'],SPEC['tube_y_lanes']):
            level=(y/radii[1])**2+(z/radii[2])**2
            if level<SPEC['layout_ellipsoid_norm_limit']:
                rows.append(dict(y_model=float(y),z_model=float(center[2]+z),
                    half_length_model=float(radii[0]*np.sqrt(1-level)+SPEC['tube_end_extension_model'])))
    return rows


def identity():
    return hashlib.sha256(json.dumps(dict(spec=SPEC,tube_layout=tube_layout()),sort_keys=True,
                                    separators=(',',':'),allow_nan=False).encode()).hexdigest()


def dense_centerline(index):
    layout=tube_layout()[index];u=np.linspace(0.,1.,SPEC['tube_samples']);phase=SPEC['tube_wave_phase_increment']*index
    x=(2*u-1)*layout['half_length_model']
    # Every tendril crosses the central vertical hub. Endpoint lane coordinates
    # fan out to fill the ellipsoidal envelope, rather than leaving floating rods.
    y=layout['y_model']*(2*u-1)**2+SPEC['tube_wave_y_model']*np.sin(4*np.pi*u)*np.sin(np.pi*u)
    z=layout['z_model']+SPEC['tube_wave_z_model']*np.sin(3*np.pi*u+phase)*np.sin(np.pi*u)
    return np.column_stack((x,y,z))


def components():
    sphere=trimesh.creation.icosphere(subdivisions=SPEC['ellipsoid_subdivisions'])
    center=np.asarray(SPEC['ellipsoid_center_model'])
    sphere.vertices=sphere.vertices*np.asarray(SPEC['ellipsoid_radii_model'])+center
    inner=sphere.copy();inner.vertices=center+SPEC['inner_scale']*(inner.vertices-center)
    hub=trimesh.creation.cylinder(radius=SPEC['hub_radius_model'],height=np.diff(SPEC['hub_interval_model'])[0],sections=24)
    hub.apply_translation([0,0,np.mean(SPEC['hub_interval_model'])])
    profile=[(0.,-.16)]+[(.025*np.cos(a),-.135+.025*np.sin(a)) for a in np.linspace(-np.pi/2,0,7)[1:]]
    r,z=SPEC['grip_radius_model'],SPEC['grip_center_z_model'];start=-np.arccos(SPEC['shaft_radius_model']/r)
    profile += [(SPEC['shaft_radius_model'],-.03),(SPEC['shaft_radius_model'],z+r*np.sin(start))]
    profile += [(r*np.cos(a),z+r*np.sin(a)) for a in np.linspace(start,np.pi/2,17)[1:-1]]+[(0.,z+r)]
    pin=trimesh.creation.revolve(np.asarray(profile),sections=32)
    if pin.volume<0:pin.invert()
    result={'outer-ellipsoid':sphere,'inner-ellipsoid-candidate':inner,'hub-fixed':hub,'pin-fixed':pin}
    for i in range(len(tube_layout())):
        result[f'tube-{i:03d}-fixed']=tube_mesh(dense_centerline(i),SPEC['tube_radius_model'],SPEC['tube_sections'])
    return result


def preview(directory,parts):
    import pybullet as p
    from PIL import Image,ImageDraw,ImageFont
    out=Image.new('RGB',(1600,1080),'#f2f4f7');draw=ImageDraw.Draw(out)
    font=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',25)
    small=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',20)
    client=p.connect(p.DIRECT)
    try:
        ids=[]
        for name,mesh in parts.items():
            if 'ellipsoid' in name:continue
            shape=p.createVisualShape(p.GEOM_MESH,fileName=str(directory/(name+'-model.obj')),
                rgbaColor=[.24,.28,.31,1] if 'tube' in name else [.62,.65,.69,1],physicsClientId=client)
            ids.append(p.createMultiBody(0,baseVisualShapeIndex=shape,baseCollisionShapeIndex=-1,physicsClientId=client))
        cameras=[([0,-5,1.1],[0,0,1.]),([5,0,1.1],[0,0,1.]),([3.7,-5,3.3],[0,0,1.])]
        projection=p.computeProjectionMatrixFOV(32,800/440,.01,20)
        labels=['正面：密集横向波浪细管','侧面：前后分层与中心连接','透视：握柄、管束与有限半径底尖']
        for i,(eye,target) in enumerate(cameras):
            pixels=p.getCameraImage(800,440,p.computeViewMatrix(eye,target,[0,0,1]),projection,
                renderer=p.ER_TINY_RENDERER,lightDirection=[-2,-3,5],physicsClientId=client)[2]
            panel=Image.fromarray(np.asarray(pixels,np.uint8).reshape(440,800,4)[:,:,:3])
            x=(i%2)*800;y=55+(i//2)*480
            out.paste(panel,(x,y));draw.text((x+18,y-35),labels[i],font=font,fill='#17364c')
    finally:p.disconnect(client)
    # A genuine fixed x=0 input section, not an optimized cut or invented cavity.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse,Rectangle,Circle
    fig,ax=plt.subplots(figsize=(8,4.4),dpi=100)
    center=np.asarray(SPEC['ellipsoid_center_model']);r=np.asarray(SPEC['ellipsoid_radii_model'])
    ax.add_patch(Ellipse((0,center[2]),2*r[1],2*r[2],facecolor='#df86b2',alpha=.3,edgecolor='#b64c83'))
    ax.add_patch(Ellipse((0,center[2]),2*r[1]*SPEC['inner_scale'],2*r[2]*SPEC['inner_scale'],
        facecolor='#df86b2',alpha=.3,edgecolor='#b64c83',linestyle='--'))
    ax.add_patch(Rectangle((-SPEC['hub_radius_model'],SPEC['hub_interval_model'][0]),2*SPEC['hub_radius_model'],
        np.diff(SPEC['hub_interval_model'])[0],facecolor='#60646b',alpha=.5))
    for i in range(len(tube_layout())):
        point=dense_centerline(i)[SPEC['tube_samples']//2]
        ax.add_patch(Circle((point[1],point[2]),SPEC['tube_radius_model'],facecolor='#272e35',alpha=.5))
    ax.set(xlim=(-1.,1.),ylim=(0.,2.1),xlabel='y [model]',ylabel='z [model]');ax.set_aspect('equal')
    ax.set_title('x=0 input section: pink candidate, gray fixed tendrils / hub')
    fig.tight_layout();fig.canvas.draw();pixels=np.asarray(fig.canvas.buffer_rgba())[:,:,:3];plt.close(fig)
    out.paste(Image.fromarray(pixels).resize((800,440)),(800,535))
    draw.text((818,500),'中心剖面：椭球候选被固定管和轴杆穿过',font=font,fill='#17364c')
    draw.text((18,1000),f'{len(tube_layout())} 根细管；只依据未优化参考图重建，隐藏几何为预登记假设。',font=small,fill='#17364c')
    draw.text((18,1035),'当前是输入部件预览；精确连接、材料分区、优化和物理仿真均待独立检查。',font=small,fill='#17364c')
    out.save(directory/'input-preview.png')


def prepare(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    parts=components();entries={}
    for name,mesh in parts.items():
        audit=mesh_audit(mesh)
        if not audit['closed'] or not audit['winding_consistent'] or audit['components']!=1 or mesh.volume<=0:
            raise ValueError('INVALID_CLOSED_SOURCE_PART:'+name)
        write_obj(directory/(name+'-model.obj'),mesh)
        si=mesh.copy();si.vertices*=SPEC['model_unit_m'];write_obj(directory/(name+'-si.obj'),si)
        entries[name]=dict(audit=audit,model_obj=name+'-model.obj',si_obj=name+'-si.obj',
            sha256={s:hashlib.sha256((directory/(name+'-'+s+'.obj')).read_bytes()).hexdigest() for s in ('model','si')})
    reference=hashlib.sha256(REFERENCE.read_bytes()).hexdigest() if REFERENCE.exists() else None
    fixed=['hub-fixed','pin-fixed']+[name for name in parts if name.startswith('tube-')]
    manifest=dict(schema='dense-noodles-source-parts-v2',spec=SPEC,spec_sha256=identity(),tube_layout=tube_layout(),
        reference_image_sha256=reference,parts=entries,
        construction=dict(fixed_parts=fixed,outer_operation='UNION(outer-ellipsoid, ALL fixed_parts)',
            design_operation='DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',
            fixed_material_operation='DIFFERENCE(outer, design)'),
        connections='EVERY_TENDRIL_CROSSES_CENTRAL_HUB; EXACT_UNION_NOT_YET_CHECKED',
        status='SOURCE_PARTS_PREPARED_PENDING_EXACT_ADMISSION',optimization_invocations=0,simulation_invocations=0,
        unknown_author_parameters=SPEC['hidden_geometry_assumptions'],
        creator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        tube_mesh_helper_sha256=hashlib.sha256(Path(__file__).with_name('noodles_input.py').read_bytes()).hexdigest())
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    preview(directory,parts)
    return manifest


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',type=Path,required=True)
    prepare(parser.parse_args().prepare)
