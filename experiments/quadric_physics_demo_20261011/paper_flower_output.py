"""OBJ derivation from an actual saved Flower solve; never invents its cut.

Preparation and exact-region audit are separate scheduled jobs because the
Shapely runtime and native CPython 3.13 extension use different interpreters.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import sys
from fractions import Fraction as F

CIRCLE_SEGMENTS=4096


def exact_zero_area_faces(mesh):
    witnesses=[]
    for i,face in enumerate(mesh.faces):
        a,b,c=[[F(float(x)) for x in mesh.vertices[int(index)]] for index in face]
        u=[x-y for x,y in zip(b,a)];v=[x-y for x,y in zip(c,a)]
        cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        if all(x==0 for x in cross):
            witnesses.append(dict(face_index=i,vertex_indices=[int(x) for x in face],
                                  vertices=mesh.vertices[face].tolist()))
    return witnesses


def write_roundtrip_obj(mesh,path):
    """17 significant digits preserve every intended finite binary64 vertex.

    Trimesh's digits=17 means decimal *places*, which loses significant bits
    in these millimeter-scale files and can collapse a cap triangle.
    """
    with Path(path).open('x') as stream:
        stream.write('# binary64 round-trip coordinates; meters\n')
        for point in mesh.vertices:
            if any(not math.isfinite(float(x)) for x in point):raise ValueError('NONFINITE_OBJ_VERTEX')
            stream.write('v '+' '.join(format(float(x),'.17g') for x in point)+'\n')
        for face in mesh.faces:
            stream.write('f '+' '.join(str(int(x)+1) for x in face)+'\n')


def constrained_boundary_prism(polygon,z_interval,meters_per_unit):
    """Retriangulate planar caps without moving or adding boundary vertices."""
    import numpy as np
    import trimesh
    from shapely.geometry import Polygon
    from shapely import constrained_delaunay_triangles
    from paper_flower_fixture import boundary_loops
    rings=boundary_loops(polygon)
    xy=np.asarray([point for ring in rings for point in ring],float)
    si=xy*meters_per_unit
    lookup={tuple(point):i for i,point in enumerate(si)}
    if len(lookup)!=len(si):raise ValueError('DISTINCT_BOUNDARY_VERTICES_COLLIDE_IN_SI')
    si_rings=[[[x*meters_per_unit,y*meters_per_unit] for x,y in ring] for ring in rings]
    geometry=Polygon(si_rings[0],si_rings[1:])
    triangles=constrained_delaunay_triangles(geometry)
    caps=[]
    for triangle in triangles.geoms:
        points=list(triangle.exterior.coords)[:-1]
        if len(points)!=3:raise ValueError('NONTRIANGULAR_CDT_FACE')
        try:face=[lookup[tuple(point)] for point in points]
        except KeyError as error:raise ValueError('CDT_INSERTED_OR_MOVED_BOUNDARY_VERTEX') from error
        a,b,c=[[F(float(value)) for value in si[index]] for index in face]
        cross=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
        if cross==0:raise ValueError('CDT_EXACT_ZERO_AREA')
        if cross<0:face=face[::-1]
        caps.append(face)
    n=len(xy);low,high=z_interval
    vertices=np.vstack([np.column_stack([xy,np.full(n,low)]),np.column_stack([xy,np.full(n,high)])])
    faces=[list(reversed(face)) for face in caps]+[[i+n for i in face] for face in caps]
    offset=0
    for ring in rings:
        for j in range(len(ring)):
            a=offset+j;b=offset+(j+1)%len(ring)
            faces.extend([[a,b,b+n],[a,b+n,a+n]])
        offset+=len(ring)
    mesh=trimesh.Trimesh(vertices=vertices,faces=faces,process=False)
    if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume<=0:
        raise ValueError('INVALID_CONSTRAINED_CAP_PRISM')
    return mesh


def profile_geometry(loops):
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    from paper_flower_fixture import exact_profile_moments
    # Frozen boundary_loops writes each outer ring immediately before its
    # holes. Keep that grouping to preserve islands inside other holes.
    parts=[];shell=None;holes=[]
    for loop in loops:
        if exact_profile_moments([loop])[0]>0:
            if shell is not None:parts.append(Polygon(shell,holes))
            shell=loop;holes=[]
        else:
            if shell is None:raise ValueError('HOLE_WITHOUT_PARENT_OUTER')
            holes.append(loop)
    if shell is not None:parts.append(Polygon(shell,holes))
    return unary_union(parts)


def shared_bearing_profiles(noded_loops,bearing_index,seat_loops_si):
    """Reuse the noded bore on both mating sides of the derived assembly."""
    from shapely.geometry import Polygon
    bore=Polygon(noded_loops[bearing_index])
    seat=Polygon(noded_loops[bearing_index],seat_loops_si[1:])
    if not bore.is_valid or not seat.is_valid or seat.area<=0:
        raise ValueError('INVALID_SHARED_BEARING_PROFILE')
    return bore,seat


def raw_mesh_moments(mesh,density=1.):
    import numpy as np
    p=mesh.mass_properties
    mass=float(p['mass']);c=np.asarray(p['center_mass'])
    central=np.asarray(p['inertia'])
    origin=central+mass*(np.dot(c,c)*np.eye(3)-np.outer(c,c))
    q=.5*np.trace(origin)*np.eye(3)-origin
    return density*np.array([mass,*(mass*c),q[0,0],q[1,1],q[2,2],q[1,2],q[0,2],q[0,1]])


def prepare(folder,application,out,*,shared_arrangement=False):
    import numpy as np
    import trimesh
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    from paper_flower_fixture import polygons,boundary_loops,exact_profile_moments
    from solid_view import extrude
    record=json.loads((folder/'input.json').read_text())
    manifest=json.loads((folder/'flower-cases.json').read_text())
    problem=manifest['cases'][0]['planar_problem']
    data=json.loads(application.read_text())
    case=next(c for c in data['cases'] if c['id']==manifest['cases'][0]['id'])
    if case['planar_problem']!=problem:raise ValueError('SAVED_REQUEST_MISMATCH')
    attempt=next((a for a in case['attempts'] if a['seed']==0 and a['status']=='CONVERGED'),None)
    if attempt is None:raise ValueError('PREREGISTERED_SEED_0_NOT_CONVERGED')
    solved=attempt['solver_result']
    if solved['accepted_steps']<=0:raise ValueError('NO_ACTUAL_OPTIMIZATION_UPDATE')
    # Reconstruct each signed planar region using the frozen loops.
    design=profile_geometry(problem['design_region_loops'])
    center=solved['center'];radius=solved['radius']
    if not math.isfinite(radius) or radius<=0:raise ValueError('INVALID_SAVED_RADIUS')
    cutter=Polygon([(center[0]+radius*math.cos(2*math.pi*i/CIRCLE_SEGMENTS),
                     center[1]+radius*math.sin(2*math.pi*i/CIRCLE_SEGMENTS)) for i in range(CIRCLE_SEGMENTS)])
    # The same cut polygon supplies both sides, rather than independent arc
    # tessellations whose vertices or chords can disagree.
    original=Polygon(record['retained_original_loops'][0],record['retained_original_loops'][1:])
    authority_span=float(max(design.bounds[2]-design.bounds[0],design.bounds[3]-design.bounds[1]))
    shared=None
    if shared_arrangement:
        from paper_flower_shared_partition import shared_partition,circle_chord_bound
        length=record['meters_per_model_unit'];ls=F(float(length))
        shared=shared_partition(original,design,cutter,length)
        cp=[tuple(F(float(x))*ls for x in p) for p in list(cutter.exterior.coords)[:-1]]
        cc=tuple(F(float(x))*ls for x in center);rr=F(float(radius))*ls
        radial=circle_chord_bound(cp,cc,rr)
        rounding=F(shared['certificate']['maximum_coordinate_rounding_l1_m'])
        bound=radial+rounding
        tube=F(shared['certificate']['design_boundary_tube_area_bound_m2'])
        # Circle annulus bounds both directions, including complete arcs;
        # the interface tube bounds the noded design/rim perturbation.
        area_bound=4*F(355,113)*rr*bound+tube
        height=F(float(problem['extrusion_interval'][1]))-F(float(problem['extrusion_interval'][0]))
        zerror=max(abs(F(float(F(float(z))*ls))-F(float(z))*ls) for z in problem['extrusion_interval'])
        bboxarea=F(float(original.bounds[2]-original.bounds[0]))*F(float(original.bounds[3]-original.bounds[1]))*ls**2
        volume_bound=area_bound*(height*ls+2*zerror)+2*zerror*bboxarea
        input_span=F(float(max(original.bounds[2]-original.bounds[0],original.bounds[3]-original.bounds[1])))*ls
        if bound>input_span*F(1,10**6) or volume_bound>input_span**3*F(1,10**6):
            raise ValueError('BOUNDED_APPROXIMATION_BUDGET_EXCEEDED')
        shared['certificate'].update(circle_center_model=center,circle_radius_model=radius,
            circle_segments=CIRCLE_SEGMENTS,input_span_m=str(input_span),
            circle_radial_bound_m=str(radial),geometric_deviation_bound_m=str(bound),
            symmetric_difference_volume_bound_m3=str(volume_bound),
            maximum_geometric_budget_m=str(input_span*F(1,10**6)),
            maximum_volume_budget_m3=str(input_span**3*F(1,10**6)),
            authority_scope='UNCHANGED_FROZEN_SOURCE_LOOPS_AND_SAVED_ANALYTIC_CIRCLE',
            derived_scope='SINGLE_ROUNDED_SHARED_PLANAR_ARRANGEMENT')
        design=shared['design'];original=shared['original']
        retained=shared['retained'];cavity=shared['cavity']
    else:
        retained=design.difference(cutter);cavity=design.intersection(cutter)
    if not retained.is_valid or not cavity.is_valid:raise ValueError('INVALID_POLYGONIZED_CUT')
    if cavity.area<=0:raise ValueError('NO_ACTUAL_MATERIAL_CHANGE')
    out.mkdir(parents=True,exist_ok=False)
    def mesh_for(geometry,z):
        if geometry.is_empty:raise ValueError('EMPTY_DERIVED_REGION')
        if shared is not None:
            z=[float(F(float(v))*F(float(record['meters_per_model_unit']))) for v in z]
        return trimesh.util.concatenate([constrained_boundary_prism(p,z,1. if shared is not None else record['meters_per_model_unit']) for p in polygons(geometry)])
    def to_si(geometry):
        if shared is None:return geometry
        length=F(float(record['meters_per_model_unit']))
        return profile_geometry([[[float(F(float(x))*length),float(F(float(y))*length)] for x,y in loop] for loop in boundary_loops(geometry)])
    def model_raw(mesh,density):
        raw=raw_mesh_moments(mesh,density)
        if shared is not None:
            length=record['meters_per_model_unit'];raw[0]/=length**3;raw[1:4]/=length**4;raw[4:]/=length**5
        return raw
    meshes={'design':mesh_for(design,problem['extrusion_interval']),
            'retained':mesh_for(retained,problem['extrusion_interval']),
            'cavity':mesh_for(cavity,problem['extrusion_interval'])}
    meshes['original-flower']=mesh_for(original,[-.08,.08])
    meshes['original-middle']=mesh_for(original,[-.06,.06])
    shared_seat=None
    if shared is not None and shared['certificate'].get('source_and_design_share_noded_boundary'):
        seat_layer=next(x for x in record['layers'] if x['id']=='bearing-seat')
        length=F(float(record['meters_per_model_unit']))
        seat_si=[[[float(F(float(x))*length),float(F(float(y))*length)] for x,y in loop] for loop in seat_layer['loops']]
        bore,shared_seat=shared_bearing_profiles(shared['certificate']['noded_design_loops_si'],record['bearing_loop_index'],seat_si)
    else:
        bore=to_si(Polygon(record['retained_original_loops'][record['bearing_loop_index']]))
    meshes['bearing-bore-region']=mesh_for(bore,[-.065,.065])
    for layer in record['layers']:
        # Rebuild from the frozen model-coordinate loops, never an SI export
        # round trip that silently perturbs the authoritative vertices.
        if shared is not None and layer['id']=='fixed-middle-boundary-rim':
            geometry=shared['rim']
        elif shared is not None and shared['certificate'].get('source_and_design_share_noded_boundary') and layer['id'] in ('top','bottom'):
            geometry=shared['original']
        elif layer['id']=='bearing-seat' and shared_seat is not None:
            geometry=shared_seat
        else:
            geometry=to_si(profile_geometry(layer['loops']))
        meshes[layer['id']]=mesh_for(geometry,layer['z_interval'])
    densitymap={layer['id']:layer['volume_density'] for layer in record['layers']}
    optimized=model_raw(meshes['retained'],problem['extrusion_volume_density'])
    original_raw=model_raw(meshes['design'],problem['extrusion_volume_density'])
    for name,density in densitymap.items():
        contribution=model_raw(meshes[name],density)
        optimized+=contribution;original_raw+=contribution
    from flower_grip_preview import assembly_spec
    from paper_flower_fixture import circle
    clamp_densitymap={}
    clamp_raw=np.zeros(10)
    for name,spec in assembly_spec().items():
        if spec['motion']!='FIXED_CLAMP':continue
        profile=circle(spec['radius'])
        if 'inner_radius' in spec:profile=profile.difference(circle(spec['inner_radius']))
        meshes[name]=mesh_for(to_si(profile),spec['z'])
        density=6.5 if name in ('shaft','inner_race') else 1.
        clamp_densitymap[name]=density
        clamp_raw+=model_raw(meshes[name],density)
    expected=np.asarray(solved['extruded_complete_material_raw_moments'])
    span=authority_span
    scale=original_raw[0]*np.array([1.]+[span]*3+[span**2]*6)
    error=np.abs(optimized-expected)/scale
    files={}
    for name,mesh in meshes.items():
        if shared is None:mesh.apply_scale(record['meters_per_model_unit'])
        before=exact_zero_area_faces(mesh)
        if before:raise ValueError('PREEXPORT_EXACT_ZERO_AREA: '+name+' '+str(before[:3]))
        path=out/(name+'.obj');write_roundtrip_obj(mesh,path)
        loaded=trimesh.load(path,force='mesh',process=False)
        after=exact_zero_area_faces(loaded)
        if after:raise ValueError('RELOAD_EXACT_ZERO_AREA: '+name+' '+str(after[:3]))
        if not np.array_equal(mesh.vertices,loaded.vertices):raise ValueError('OBJ_BINARY64_ROUNDTRIP_CHANGED: '+name)
        files[name]=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         watertight=bool(loaded.is_watertight),winding_consistent=bool(loaded.is_winding_consistent),
                         positive_volume=bool(loaded.volume>0),volume_m3=float(loaded.volume),
                         exact_zero_area_face_count=0,binary64_roundtrip_vertices_unchanged=True)
    def si_properties(raw):
        v=raw.copy();length=record['meters_per_model_unit'];density=record['pla_density_kg_m3']
        v[0]*=density*length**3;v[1:4]*=density*length**4;v[4:]*=density*length**5
        m=v[0];c=v[1:4]/m;q=np.array([[v[4],v[9],v[8]],[v[9],v[5],v[7]],[v[8],v[7],v[6]]])
        inertia=np.trace(q)*np.eye(3)-q-m*(np.dot(c,c)*np.eye(3)-np.outer(c,c))
        return dict(mass_kg=float(m),center_m=c.tolist(),inertia_kg_m2=inertia.tolist())
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'python'))
    from make_it_stand.verification.planar_material import integrate_planar_entity
    authority=np.sum([np.asarray(integrate_planar_entity(e,origin=[0,0],length_scale=span)['raw_mass_moments'])
                      for e in solved['extruded_material_entities']],axis=0)
    authority_error=np.abs(authority-expected)/scale
    result=dict(schema='paper-flower-derived-output-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',
        input_sha256=hashlib.sha256((folder/'input.json').read_bytes()).hexdigest(),
        application_sha256=hashlib.sha256(application.read_bytes()).hexdigest(),selected_seed=0,
        circle_center=center,circle_radius=radius,circle_segments=CIRCLE_SEGMENTS,
        maximum_model_sagitta=radius*(1-math.cos(math.pi/CIRCLE_SEGMENTS)),
        circle_representation='INSCRIBED_POLYGON_DISPLAY_AND_OBJ_DERIVATIVE_NOT_AUTHORITY_ARC',
        obj_coordinate_encoding='17_SIGNIFICANT_DIGIT_BINARY64_ROUNDTRIP',
        cap_triangulation='GEOS_CONSTRAINED_DELAUNAY_NO_BOUNDARY_VERTEX_MOVEMENT_OR_INSERTION',
        shared_cut_polygon=True,planar_partition_symmetric_difference_area=float(design.symmetric_difference(retained.union(cavity)).area),
        planar_overlap_area=float(retained.intersection(cavity).area),
        independent_triangle_moment_error=error.tolist(),moment_gate_pass=bool(np.max(error)<=1e-6),
        normalization_initial_model_mass=float(original_raw[0]),normalization_input_span=span,
        authority_independent_moment_error=authority_error.tolist(),
        original=si_properties(original_raw),optimized=si_properties(authority),
        derived_optimized=si_properties(optimized),
        clamp=si_properties(clamp_raw),clamp_layers=clamp_densitymap,
        clamp_density_scope='ENGINEERING_ASSUMPTION_PLA_CAPS_STEEL_SHAFT_AND_INNER_RACE',
        files=files,rotor_layers=densitymap,excluded_from_rotor=record['excluded_from_rotor'],
        full_material_gate='PENDING_EXACT_REGION_AUDIT',simulation_admission=False)
    if shared is not None:
        attachment=out/'shared-arrangement.json'
        attachment.write_text(json.dumps(shared['certificate'],indent=2,allow_nan=False)+'\n')
        result.update(derived_geometry_policy='BOUNDED_SHARED_RATIONAL_ARRANGEMENT_V1',
                      shared_arrangement_attachment=attachment.name,
                      shared_arrangement_sha256=hashlib.sha256(attachment.read_bytes()).hexdigest())
    (out/'output-preparation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['moment_gate_pass','independent_triangle_moment_error','maximum_model_sagitta','full_material_gate']}))


def _exact_audit_impl(out,report_directory,event):
    import numpy as np
    import trimesh
    import _make_it_stand as native
    record=json.loads((out/'output-preparation.json').read_text())
    regions={}
    for info in record['files'].values():
        if hashlib.sha256((out/info['file']).read_bytes()).hexdigest()!=info['sha256']:
            raise ValueError('DERIVED_OBJ_HASH_MISMATCH: '+info['file'])
    for name in ['design','retained','cavity','original-middle','original-flower','bearing-bore-region',*record['rotor_layers']]:
        event(dict(stage='authenticate_operand',name=name,file=str(out/(name+'.obj'))))
        mesh=trimesh.load(out/(name+'.obj'),force='mesh',process=False)
        regions[name]=native.make_exact_material_region(np.asarray(mesh.vertices),np.asarray(mesh.faces,dtype=np.int32))
        event(dict(stage='operand_authenticated',name=name))
    event(dict(stage='editable_partition'))
    partition=dict(native.check_exact_material_partition(regions['design'],[regions['retained'],regions['cavity']]))
    event(dict(stage='middle_partition'))
    middle_partition=dict(native.check_exact_material_partition(regions['original-middle'],[regions['design'],regions['fixed-middle-boundary-rim']]))
    event(dict(stage='original_partition'))
    original_partition=dict(native.check_exact_material_partition(regions['original-flower'],[regions['bottom'],regions['top'],regions['fixed-middle-boundary-rim'],regions['design']]))
    event(dict(stage='bearing_containment'))
    bore_subsets={name:dict(native.check_exact_material_subset(regions[name],regions['bearing-bore-region'])) for name in ['bearing-seat','outer-race']}
    pairs=[]
    names=['retained',*record['rotor_layers']]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            event(dict(stage='material_pair',a=a,b=b))
            check=dict(native.relate_exact_material_regions(regions[a],regions[b]))
            pairs.append(dict(a=a,b=b,relation=check))
    # This historical entry point checks represented sets only. Producer
    # moment flags are not independent evidence and cannot admit simulation.
    passed=partition['completed'] and partition['complete_partition'] and partition['no_positive_volume_overlap']
    passed=passed and all(p['completed'] and p['complete_partition'] and p['no_positive_volume_overlap'] for p in [middle_partition,original_partition])
    passed=passed and all(p['completed'] and p['a_subset_b'] for p in bore_subsets.values())
    passed=passed and all(row['relation']['completed'] and not row['relation']['positive_volume_overlap'] for row in pairs)
    result=dict(schema='paper-flower-derived-exact-audit-v1',preparation_sha256=hashlib.sha256((out/'output-preparation.json').read_bytes()).hexdigest(),
                native_module=str(native.__file__),native_module_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                partition=partition,middle_partition=middle_partition,original_partition=original_partition,
                bore_subsets=bore_subsets,pairs=pairs,full_material_gate='UNKNOWN' if passed else 'BLOCKED',
                geometry_partition_checks_passed=bool(passed),
                admission_policy='LEGACY_GEOMETRY_ONLY_NON_ADMITTING_V2',
                required_final_checker='paper_flower_saved_check.py',
                simulation_admission=False,scope='EXACT_REPRESENTED_POLYHEDRAL_OBJ_SETS_ONLY_NO_INDEPENDENT_MOMENT_OR_SOURCE_ADMISSION')
    report_directory.mkdir(parents=True,exist_ok=True)
    with (report_directory/'exact-audit.json').open('x') as f:f.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['full_material_gate','simulation_admission']}))
    event(dict(stage='audit_completed',full_material_gate=result['full_material_gate']))


def exact_audit(out,report_directory):
    report_directory.mkdir(parents=True,exist_ok=True)
    last={}
    with (report_directory/'audit-events.jsonl').open('x') as journal:
        def event(row):
            last.clear();last.update(row)
            journal.write(json.dumps(row,allow_nan=False)+'\n');journal.flush();os.fsync(journal.fileno())
        event(dict(stage='audit_started',output=str(out)))
        try:
            _exact_audit_impl(out,report_directory,event)
        except Exception as error:
            failure=dict(schema='paper-flower-derived-audit-failure-v1',
                         full_material_gate='BLOCKED',simulation_admission=False,
                         last_completed_or_started_event=dict(last),exception_type=type(error).__name__,
                         exception_message=str(error),preparation_sha256=hashlib.sha256((out/'output-preparation.json').read_bytes()).hexdigest())
            native=sys.modules.get('_make_it_stand')
            if native is not None:
                failure['native_module']=str(native.__file__)
                failure['native_module_sha256']=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()
            with (report_directory/'audit-failure.json').open('x') as attachment:
                attachment.write(json.dumps(failure,indent=2,allow_nan=False)+'\n');attachment.flush();os.fsync(attachment.fileno())
            event(dict(stage='audit_failed',exception_type=type(error).__name__,exception_message=str(error)))
            raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--folder',type=Path);parser.add_argument('--application',type=Path);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--exact-audit',action='store_true');parser.add_argument('--audit-report-directory',type=Path);parser.add_argument('--shared-arrangement',action='store_true');args=parser.parse_args()
    if args.exact_audit:exact_audit(args.output,args.audit_report_directory or args.output)
    else:prepare(args.folder,args.application,args.output,shared_arrangement=args.shared_arrangement)
