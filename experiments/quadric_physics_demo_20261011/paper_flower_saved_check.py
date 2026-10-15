"""Independent saved-object gate; never trusts output preparation pass flags.

Runs only as a separately scheduled job. Input, application, approved clamp
reference and loaded native module must be bound by the scheduling receipt.
"""
from pathlib import Path
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.metadata
import json
import math
import os
import sys
import numpy as np
import trimesh


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition,message):
    if not condition:raise ValueError(message)


def verify_native_identity(native,expected_sha):
    require(len(expected_sha)==64 and sha(native.__file__)==expected_sha,'NATIVE_BINARY_IDENTITY_MISMATCH')
    return dict(path=str(native.__file__),sha256=expected_sha)


def validate_run_identity(case,problem,prep):
    require(case['planar_problem']==problem,'FROZEN_REQUEST_MISMATCH')
    require(prep['selected_seed']==0,'PREDECLARED_SEED_MISMATCH')
    attempt=next((a for a in case['attempts'] if a['seed']==0 and a['status']=='CONVERGED'),None)
    require(attempt is not None,'NO_PREREGISTERED_SUCCESS')
    solved=attempt['solver_result']
    require(solved['accepted_steps']>0,'NO_ACTUAL_ACCEPTED_UPDATE')
    require(solved['center']==prep['circle_center'] and solved['radius']==prep['circle_radius'],'SAVED_CIRCLE_IDENTITY_MISMATCH')
    require(solved['radius']>0 and math.isfinite(solved['radius']),'INVALID_RADIUS')
    return solved


def integrate_saved_mesh(mesh,density):
    """Signed tetrahedron simplex integrals; independent of mesh.mass_properties."""
    require(math.isfinite(density) and density>0,'INVALID_DENSITY')
    rows=[[] for _ in range(10)]
    for face in mesh.faces:
        a,b,c=np.asarray(mesh.vertices[face],float)
        volume=float(np.dot(a,np.cross(b,c)))/6.
        s=a+b+c
        q=volume*(np.outer(s,s)+np.outer(a,a)+np.outer(b,b)+np.outer(c,c))/20.
        first=volume*s/4.
        terms=[volume,*first,q[0,0],q[1,1],q[2,2],q[1,2],q[0,2],q[0,1]]
        for row,value in zip(rows,terms):row.append(float(value))
    return density*np.asarray([math.fsum(row) for row in rows])


def expected_wall_triangles(loops,z_interval,length):
    low,high=[float(z)*length for z in z_interval]
    result=Counter()
    for loop in loops:
        for a,b in zip(loop,loop[1:]+loop[:1]):
            aa=(float(a[0])*length,float(a[1])*length,low)
            bb=(float(b[0])*length,float(b[1])*length,low)
            at=(*aa[:2],high);bt=(*bb[:2],high)
            result[tuple(sorted((aa,bb,bt)))]+=1
            result[tuple(sorted((aa,bt,at)))]+=1
    return result


def check_frozen_prism_boundary(mesh,loops,z_interval,length):
    low,high=[float(z)*length for z in z_interval]
    expected={(float(p[0])*length,float(p[1])*length,z) for loop in loops for p in loop for z in (low,high)}
    require({tuple(v) for v in mesh.vertices}==expected,'FROZEN_PRISM_VERTICES_CHANGED')
    walls=Counter()
    for face in mesh.faces:
        points=[tuple(v) for v in mesh.vertices[face]]
        if not(all(p[2]==low for p in points) or all(p[2]==high for p in points)):
            walls[tuple(sorted(points))]+=1
    require(walls==expected_wall_triangles(loops,z_interval,length),'FROZEN_PRISM_BOUNDARY_CHANGED')


def properties(raw):
    mass=raw[0];require(mass>0,'NONPOSITIVE_MASS')
    c=raw[1:4]/mass
    q=np.array([[raw[4],raw[9],raw[8]],[raw[9],raw[5],raw[7]],[raw[8],raw[7],raw[6]]])
    inertia=np.trace(q)*np.eye(3)-q-mass*(np.dot(c,c)*np.eye(3)-np.outer(c,c))
    require(np.linalg.eigvalsh(inertia).min()>0,'NONPOSITIVE_INERTIA')
    return dict(mass_kg=float(mass),center_m=c.tolist(),inertia_kg_m2=inertia.tolist())


def model_to_si(raw,length,density):
    return np.asarray(raw,float)*density*np.array([length**3]+[length**4]*3+[length**5]*6)


def authenticate_source_loop_order(serialized,frozen_loops):
    """Authenticate directed source edges while preserving local ring indices."""
    from fractions import Fraction as F
    ordered=[[tuple(F(x) for x in p) for p in loop] for loop in serialized]
    frozen=[[tuple(F(float(x)) for x in p) for p in loop] for loop in frozen_loops]
    edges=lambda loops:Counter((a,b) for loop in loops for a,b in zip(loop,loop[1:]+loop[:1]))
    require(all(len(loop)>=3 for loop in ordered) and edges(ordered)==edges(frozen),'SOURCE_ORDER_BOUNDARY_CHANGED')
    return ordered


def verify_shared_certificate(certificate,frozen,problem,solved):
    """Recompute error bounds and source membership, without producer flags.

Native solids still decide the represented partition. This attachment proves
only the approved approximation relation; it cannot grant geometry admission.
    """
    from fractions import Fraction as F
    require(certificate['schema']=='shared-rational-arrangement-v2','UNSUPPORTED_SHARED_CERTIFICATE')
    require(certificate['circle_center_model']==solved['center'] and certificate['circle_radius_model']==solved['radius'],'CERTIFICATE_CIRCLE_CHANGED')
    require(certificate['circle_segments']==4096,'CERTIFICATE_RESOLUTION_CHANGED')
    length=F(float(frozen['meters_per_model_unit']))
    source=authenticate_source_loop_order(certificate['authority_original_loops_model'],frozen['retained_original_loops'])
    design=authenticate_source_loop_order(certificate['authority_design_loops_model'],problem['design_region_loops'])
    center=solved['center'];radius=solved['radius']
    circle=[tuple(F(float(x)) for x in (center[0]+radius*math.cos(2*math.pi*i/4096),center[1]+radius*math.sin(2*math.pi*i/4096))) for i in range(4096)]
    edges=[(a,b) for loop in design for a,b in zip(loop,loop[1:]+loop[:1])]
    cutedges=list(zip(circle,circle[1:]+circle[:1]))
    nodes=[{a,b} for a,b in edges];cutnodes=[{a,b} for a,b in cutedges]
    allowed={p for loop in source+design for p in loop}|set(circle)
    for row in certificate['intersections']:
        point=tuple(F(x) for x in row['point']);a,b=edges[row['design_edge']];c,d=cutedges[row['circle_edge']]
        for start,end in [(a,b),(c,d)]:
            u=tuple(x-y for x,y in zip(end,start));v=tuple(x-y for x,y in zip(point,start))
            require(u[0]*v[1]-u[1]*v[0]==0 and all(min(x,y)<=z<=max(x,y) for x,y,z in zip(start,end,point)),'INVALID_RATIONAL_INTERSECTION_SOURCE')
        nodes[row['design_edge']].add(point);allowed.add(point)
        cutnodes[row['circle_edge']].add(point)
    mapping={};rounding=F(0)
    for row in certificate['vertices']:
        p=tuple(F(x) for x in row['authority_model']);stored=tuple(row['stored_si'])
        require(p in allowed,'UNAUTHENTICATED_SHARED_VERTEX')
        require(stored==tuple(float(x*length) for x in p),'NOT_SINGLE_ROUNDED_AUTHORITY_CONVERSION')
        require(p not in mapping,'DUPLICATE_SHARED_VERTEX')
        mapping[p]=stored
        rounding=max(rounding,sum(abs(F(x)-p[i]*length) for i,x in enumerate(stored)))
    require(set(mapping)==allowed,'MISSING_SHARED_AUTHORITY_VERTEX')
    require(len(set(mapping.values()))==len(mapping),'COLLAPSED_SHARED_AUTHORITY_VERTEX')
    noded=[];offset=0
    for loop in design:
        ring=[]
        for _ in loop:
            a,b=edges[offset];axis=0 if a[0]!=b[0] else 1
            chain=sorted(nodes[offset],key=lambda p:(p[axis]-a[axis])/(b[axis]-a[axis]))
            ring.extend(mapping[p] for p in chain[:-1]);offset+=1
        noded.append(ring)
    require(noded==[[tuple(p) for p in loop] for loop in certificate['noded_design_loops_si']],'NODED_DESIGN_IDENTITY_MISMATCH')
    design_walls={(tuple(a),tuple(b)) for loop in noded for a,b in zip(loop,loop[1:]+loop[:1])}
    original_walls={(mapping[a],mapping[b]) for loop in source for a,b in zip(loop,loop[1:]+loop[:1])}
    cut_walls=set()
    for (a,b),points in zip(cutedges,cutnodes):
        axis=0 if a[0]!=b[0] else 1
        chain=sorted(points,key=lambda p:(p[axis]-a[axis])/(b[axis]-a[axis]))
        cut_walls.update((mapping[p],mapping[q]) for p,q in zip(chain,chain[1:]))
    reverse=lambda pairs:{(b,a) for a,b in pairs}
    walls=dict(design=design_walls,retained=design_walls|reverse(cut_walls),
               cavity=design_walls|cut_walls,
               **{'fixed-middle-boundary-rim':original_walls|reverse(design_walls)})
    cc=tuple(F(float(x))*length for x in center);rr=F(float(radius))*length;radial=F(0)
    for p,q in cutedges:
        a=tuple(p[i]*length-cc[i] for i in range(2));d=tuple((q[i]-p[i])*length for i in range(2));b=tuple(q[i]*length-cc[i] for i in range(2))
        require(a[0]*b[1]-a[1]*b[0]>0,'CIRCLE_ANGULAR_ORDER_FAILURE')
        dd=sum(x*x for x in d);require(dd>0,'ZERO_CIRCLE_EDGE')
        t=max(F(0),min(F(1),-sum(a[i]*d[i] for i in range(2))/dd))
        qvalues=[sum(x*x for x in a),sum(x*x for x in b),sum((a[i]+t*d[i])**2 for i in range(2))]
        radial=max(radial,max(abs(v-rr**2)/rr for v in qvalues))
    boundary=edges+[(a,b) for loop in source for a,b in zip(loop,loop[1:]+loop[:1])]
    perimeter=sum(abs(a[0]-b[0])+abs(a[1]-b[1]) for a,b in boundary)*length
    tube=2*rounding*perimeter+F(355,113)*rounding**2*len(boundary)
    bound=radial+rounding;area=4*F(355,113)*rr*bound+tube
    low,high=map(lambda z:F(float(z)),problem['extrusion_interval'])
    zerror=max(abs(F(float(z*length))-z*length) for z in (low,high))
    allpoints=[p for loop in source for p in loop]
    # Match the frozen binary64 input span contract, not a producer-selected L.
    extent=[float(max(p[i] for p in allpoints))-float(min(p[i] for p in allpoints)) for i in range(2)]
    span=F(float(max(extent)))*length;bboxarea=F(float(extent[0]))*F(float(extent[1]))*length**2
    volume=area*((high-low)*length+2*zerror)+2*zerror*bboxarea
    checks={'maximum_coordinate_rounding_l1_m':rounding,'design_boundary_tube_area_bound_m2':tube,
            'circle_radial_bound_m':radial,'geometric_deviation_bound_m':bound,
            'symmetric_difference_volume_bound_m3':volume,'input_span_m':span}
    require(all(F(certificate[k])==v for k,v in checks.items()),'APPROXIMATION_BOUND_RECOMPUTATION_MISMATCH')
    require(bound<=span*F(1,10**6) and volume<=span**3*F(1,10**6),'APPROXIMATION_BOUND_EXCEEDED')
    return dict(independently_recomputed=True,geometric_bound_m=float(bound),volume_bound_m3=float(volume),noded_design_loops_si=noded,
                authority_vertices_si=set(mapping.values()),allowed_directed_wall_edges=walls)


def check_shared_wall_sources(mesh,name,proof,z_interval,length):
    low,high=[z*length for z in z_interval]
    require(low<high,'INVALID_SHARED_WALL_HEIGHT')
    def oriented_triangle(points):
        # Cyclic changes preserve winding; reversal is deliberately distinct.
        return min(tuple(points[i:]+points[:i]) for i in range(3))
    walls={}
    for face in mesh.faces:
        points=[tuple(v) for v in mesh.vertices[face]]
        if all(p[2]==low for p in points) or all(p[2]==high for p in points):continue
        require(all(p[2] in (low,high) for p in points),'UNREGISTERED_SHARED_WALL_HEIGHT')
        xy={tuple(p[:2]) for p in points}
        require(len(xy)==2,'SHARED_WALL_NOT_VERTICAL_EXTRUSION')
        if sum(p[2]==low for p in points)==2:
            edge=next((tuple(a[:2]),tuple(b[:2])) for a,b in zip(points,points[1:]+points[:1]) if a[2]==low and b[2]==low)
        else:
            topedge=next((tuple(a[:2]),tuple(b[:2])) for a,b in zip(points,points[1:]+points[:1]) if a[2]==high and b[2]==high)
            edge=tuple(reversed(topedge))
        require(edge in proof['allowed_directed_wall_edges'][name],'SHARED_WALL_NOT_AUTHENTICATED_SOURCE_OR_CUT')
        walls.setdefault(edge,Counter())[oriented_triangle(points)]+=1
    for (a,b),actual in walls.items():
        aa=(*a,low);bb=(*b,low);at=(*a,high);bt=(*b,high)
        first=Counter([oriented_triangle([aa,bb,bt]),oriented_triangle([aa,bt,at])])
        second=Counter([oriented_triangle([aa,bb,at]),oriented_triangle([bb,bt,at])])
        require(actual==first or actual==second,'SHARED_WALL_HALF_MISSING_DUPLICATED_OR_TWISTED')


def dependency_lock(native):
    from make_it_stand.verification import planar_material
    imported={}
    for name,module in [('checker',sys.modules[__name__]),('numpy',np),('trimesh',trimesh),('planar_material',planar_material)]:
        file=Path(module.__file__)
        imported[name]=dict(path=str(file),sha256=sha(file))
    libraries={}
    for line in Path('/proc/self/maps').read_text().splitlines():
        fields=line.split()
        if len(fields)<6:continue
        path=Path(fields[-1])
        if path.is_file() and any(token in path.name for token in ['libgmp','libmpfr','libstdc++','libpython']):
            libraries[str(path)]=sha(path)
    return dict(python=sys.version,native=native,
                versions={name:importlib.metadata.version(name) for name in ['numpy','trimesh']},
                imported_sources=imported,loaded_shared_libraries=libraries)


def run_check(folder,application,derived,assembly_reference,expected_native_sha,out,preflight_only=False):
    identity=dict(expected_sha256=expected_native_sha)
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'python'))
    from make_it_stand.verification.planar_material import integrate_planar_entity
    from paper_flower_fixture import exact_extrusion_moments
    out.mkdir(parents=True,exist_ok=True)
    events=(out/'saved-check-events.jsonl').open('x')
    stage={}
    def event(row):
        stage.clear();stage.update(row)
        events.write(json.dumps(row,allow_nan=False)+'\n');events.flush();os.fsync(events.fileno())
    try:
        event(dict(stage='native_identity'))
        import _make_it_stand as native
        identity=verify_native_identity(native,expected_native_sha)
        event(dict(stage='source_identity'))
        frozen=json.loads((folder/'input.json').read_text())
        manifest=json.loads((folder/'flower-cases.json').read_text())
        index=json.loads((folder/'input-sha256.json').read_text())
        for filename in ['input.json','flower-cases.json']:
            require(index[filename]==sha(folder/filename),'FROZEN_FILE_HASH_MISMATCH: '+filename)
        prep=json.loads((derived/'output-preparation.json').read_text())
        require(prep['input_sha256']==sha(folder/'input.json'),'PREPARATION_INPUT_IDENTITY_MISMATCH')
        require(prep['application_sha256']==sha(application),'PREPARATION_APPLICATION_IDENTITY_MISMATCH')
        problem=manifest['cases'][0]['planar_problem']
        case=next(c for c in json.loads(application.read_text())['cases'] if c['id']==manifest['cases'][0]['id'])
        solved=validate_run_identity(case,problem,prep)
        require([{k:v for k,v in row.items() if k!='id'} for row in frozen['layers']]==problem['fixed_extrusion_layers'],'FIXED_LAYER_REQUEST_IDENTITY')
        assembly=json.loads(assembly_reference.read_text())
        clamp_roles=['upper_cap','lower_cap','shaft','inner_race']
        require(set(frozen['excluded_from_rotor'])==set(clamp_roles),'FIXED_CLAMP_ROLE_MISMATCH')
        for name in clamp_roles:
            require(assembly['parts'][name]['motion']=='FIXED_CLAMP','REFERENCE_CLAMP_ROLE_MISMATCH')
        density=frozen['pla_density_kg_m3'];length=frozen['meters_per_model_unit']
        require(density>0 and length>0,'INVALID_FROZEN_UNIT_OR_DENSITY')
        densities={layer['id']:density*layer['volume_density'] for layer in frozen['layers']}
        clamp_densities={name:frozen['outer_race_density_kg_m3'] if name in ('shaft','inner_race') else density for name in clamp_roles}
        require(prep['rotor_layers']=={layer['id']:layer['volume_density'] for layer in frozen['layers']},'PREPARATION_ROTOR_DENSITY_MISMATCH')
        require(prep['clamp_layers']=={name:value/density for name,value in clamp_densities.items()},'PREPARATION_CLAMP_DENSITY_MISMATCH')
        meshes={};regions={};raw={};facts={}
        names=['design','retained','cavity','original-flower','original-middle','bearing-bore-region',*densities,*clamp_roles]
        for name in names:
            event(dict(stage='load_and_authenticate_obj',name=name))
            info=prep['files'][name];path=derived/info['file']
            require(path.parent.resolve()==derived.resolve(),'OBJ_PATH_ESCAPE')
            require(sha(path)==info['sha256'],'OBJ_HASH_MISMATCH: '+name)
            mesh=trimesh.load(path,force='mesh',process=False)
            require(mesh.is_watertight and mesh.is_winding_consistent and mesh.volume>0,'INVALID_OBJ_TOPOLOGY: '+name)
            meshes[name]=mesh
            raw[name]=integrate_saved_mesh(mesh,densities.get(name,clamp_densities.get(name,density)))
        event(dict(stage='frozen_boundary_and_units'))
        approximation=None
        if prep.get('derived_geometry_policy')=='BOUNDED_SHARED_RATIONAL_ARRANGEMENT_V1':
            certificate_path=derived/prep['shared_arrangement_attachment']
            require(certificate_path.parent.resolve()==derived.resolve() and sha(certificate_path)==prep['shared_arrangement_sha256'],'SHARED_CERTIFICATE_IDENTITY_MISMATCH')
            approximation=verify_shared_certificate(json.loads(certificate_path.read_text()),frozen,problem,solved)
            check_frozen_prism_boundary(meshes['design'],approximation['noded_design_loops_si'],[z*length for z in problem['extrusion_interval']],1.)
            for name in ['design','retained','cavity',*(['fixed-middle-boundary-rim'] if 'fixed-middle-boundary-rim' in meshes else [])]:
                require(all(tuple(v[:2]) in approximation['authority_vertices_si'] for v in meshes[name].vertices),'DERIVED_VERTEX_NOT_AUTHORITY_BOUND')
                check_shared_wall_sources(meshes[name],name,approximation,problem['extrusion_interval'],length)
        else:check_frozen_prism_boundary(meshes['design'],problem['design_region_loops'],problem['extrusion_interval'],length)
        shared_full_middle=approximation is not None and frozen.get('schema')=='paper-flower-method-input-v2'
        if shared_full_middle:
            require(not any(x['id']=='fixed-middle-boundary-rim' for x in frozen['layers']),'FULL_MIDDLE_HAS_FIXED_WEB')
            # Independently authenticate directed source edges; a caller flag
            # cannot grant the noded cover policy to a different design.
            edge_set=lambda loops:{(tuple(a),tuple(b)) for loop in loops for a,b in zip(loop,loop[1:]+loop[:1])}
            require(edge_set(frozen['retained_original_loops'])==edge_set(problem['design_region_loops']),'FULL_MIDDLE_SOURCE_DESIGN_MISMATCH')
        for layer in frozen['layers']:
            if approximation is not None and layer['id']=='fixed-middle-boundary-rim':continue
            if shared_full_middle and layer['id'] in ('top','bottom'):
                check_frozen_prism_boundary(meshes[layer['id']],approximation['noded_design_loops_si'],[z*length for z in layer['z_interval']],1.)
            elif shared_full_middle and layer['id']=='bearing-seat':
                source_bore=frozen['retained_original_loops'][frozen['bearing_loop_index']]
                undirected=lambda ring:{tuple(sorted((tuple(a),tuple(b)))) for a,b in zip(ring,ring[1:]+ring[:1])}
                require(undirected(source_bore)==undirected(layer['loops'][0]),'BEARING_SEAT_SOURCE_BORE_MISMATCH')
                noded_bore=approximation['noded_design_loops_si'][frozen['bearing_loop_index']]
                inner=[[[float(F(float(x))*F(float(length))),float(F(float(y))*F(float(length)))] for x,y in loop] for loop in layer['loops'][1:]]
                expected=[list(reversed(noded_bore)),*inner]
                check_frozen_prism_boundary(meshes[layer['id']],expected,[z*length for z in layer['z_interval']],1.)
            else:check_frozen_prism_boundary(meshes[layer['id']],layer['loops'],layer['z_interval'],length)
        if shared_full_middle:
            check_frozen_prism_boundary(meshes['bearing-bore-region'],[list(reversed(approximation['noded_design_loops_si'][frozen['bearing_loop_index']]))],[-.065*length,.065*length],1.)
        for name,interval in [('original-flower',[-.08,.08]),('original-middle',[-.06,.06])]:
            if shared_full_middle:check_frozen_prism_boundary(meshes[name],approximation['noded_design_loops_si'],[z*length for z in interval],1.)
            else:check_frozen_prism_boundary(meshes[name],frozen['retained_original_loops'],interval,length)
        for name in clamp_roles:
            s=assembly['parts'][name]
            def ring(radius):return [[radius*math.cos(2*math.pi*i/96),radius*math.sin(2*math.pi*i/96)] for i in range(96)]
            loops=[ring(s['radius'])]
            if 'inner_radius' in s:loops.append(ring(s['inner_radius'])[::-1])
            check_frozen_prism_boundary(meshes[name],loops,s['z'],length)
        event(dict(stage='independent_mass_moments'))
        original_model=sum((np.asarray([float(x) for x in exact_extrusion_moments(layer)]) for layer in frozen['layers']),start=np.zeros(10))
        original_model+=np.asarray([float(x) for x in exact_extrusion_moments(dict(loops=problem['design_region_loops'],z_interval=problem['extrusion_interval'],volume_density=problem['extrusion_volume_density']))])
        original_expected=model_to_si(original_model,length,density)
        original_actual=raw['design']+sum((raw[name] for name in densities),start=np.zeros(10))
        optimized_actual=raw['retained']+sum((raw[name] for name in densities),start=np.zeros(10))
        expected_model=np.asarray(solved['extruded_complete_material_raw_moments'])
        expected=model_to_si(expected_model,length,density)
        points=np.asarray([p for loop in problem['design_region_loops'] for p in loop])
        span=float(np.ptp(points,axis=0).max())*length
        scale=original_expected[0]*np.array([1.]+[span]*3+[span**2]*6)
        errors=dict(original=np.abs(original_actual-original_expected)/scale,
                    optimized=np.abs(optimized_actual-expected)/scale)
        independent_model=sum((np.asarray(integrate_planar_entity(e,origin=[0,0],length_scale=span/length)['raw_mass_moments']) for e in solved['extruded_material_entities']),start=np.zeros(10))
        errors['authority']=np.abs(model_to_si(independent_model,length,density)-expected)/scale
        require(all(np.isfinite(row).all() and float(row.max())<=1e-6 for row in errors.values()),'SAVED_INDEPENDENT_MOMENT_GATE_FAILURE')
        if preflight_only:
            partial=dict(schema='paper-flower-saved-independent-preflight-v1',
                         sources={str(p):sha(p) for p in [folder/'input.json',folder/'flower-cases.json',folder/'input-sha256.json',application,derived/'output-preparation.json',assembly_reference]},
                         moment_errors={name:row.tolist() for name,row in errors.items()},
                         original=properties(original_actual),optimized=properties(model_to_si(independent_model,length,density)),
                         derived_optimized=properties(optimized_actual),
                         clamp=properties(sum((raw[name] for name in clamp_roles),start=np.zeros(10))),
                         dependency_lock=dependency_lock(identity),simulation_admission=False,
                         full_material_gate='UNKNOWN',geometry_check='NOT_EXECUTED_IN_PREFLIGHT',
                         scope='SAVED_SOURCE_IDENTITY_UNITS_BOUNDARY_MOMENTS_ONLY')
            with (out/'preflight-check.json').open('x') as stream:
                stream.write(json.dumps(partial,indent=2,allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
            event(dict(stage='preflight_completed',simulation_admission=False))
            print(json.dumps({k:partial[k] for k in ['full_material_gate','simulation_admission','moment_errors']}))
            return
        for name,mesh in meshes.items():
            event(dict(stage='authenticate_exact_region',name=name))
            region=native.make_exact_material_region(np.asarray(mesh.vertices),np.asarray(mesh.faces,dtype=np.int32))
            regions[name]=region;facts[name]=dict(region.facts)
        event(dict(stage='actual_partition_and_overlap'))
        middle_parts=[regions['design']]
        if 'fixed-middle-boundary-rim' in regions:middle_parts.append(regions['fixed-middle-boundary-rim'])
        partitions=dict(editable=dict(native.check_exact_material_partition(regions['design'],[regions['retained'],regions['cavity']])),
                        middle=dict(native.check_exact_material_partition(regions['original-middle'],middle_parts)),
                        original=dict(native.check_exact_material_partition(regions['original-flower'],[regions['bottom'],regions['top'],*middle_parts])))
        subsets={name:dict(native.check_exact_material_subset(regions[name],regions['bearing-bore-region'])) for name in ['outer-race','bearing-seat']}
        overlap=[]
        physical=['retained',*densities,*clamp_roles]
        for i,a in enumerate(physical):
            for b in physical[i+1:]:
                event(dict(stage='check_material_pair',a=a,b=b))
                relation=dict(native.relate_exact_material_regions(regions[a],regions[b]))
                overlap.append(dict(a=a,b=b,relation=relation))
        passed=all(p['completed'] and p['complete_partition'] and p['no_positive_volume_overlap'] for p in partitions.values())
        passed=passed and all(p['completed'] and p['a_subset_b'] for p in subsets.values())
        passed=passed and all(p['relation']['completed'] and not p['relation']['positive_volume_overlap'] for p in overlap)
        result=dict(schema='paper-flower-saved-independent-check-v1',
                    sources={str(p):sha(p) for p in [folder/'input.json',folder/'flower-cases.json',folder/'input-sha256.json',application,derived/'output-preparation.json',assembly_reference]},
                    selected_seed=0,accepted_steps=solved['accepted_steps'],circle_center=solved['center'],circle_radius=solved['radius'],
                    moment_errors={name:row.tolist() for name,row in errors.items()},
                    original=properties(original_actual),optimized=properties(model_to_si(independent_model,length,density)),
                    derived_optimized=properties(optimized_actual),clamp=properties(sum((raw[name] for name in clamp_roles),start=np.zeros(10))),
                    partitions=partitions,bearing_subsets=subsets,material_pairs=overlap,exact_region_facts=facts,
                    dependency_lock=dependency_lock(identity),simulation_admission=bool(passed),full_material_gate='PASSED' if passed else 'BLOCKED',
                    manufacturing_approved=False)
        with (out/'saved-independent-check.json').open('x') as stream:
            stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
        event(dict(stage='completed',full_material_gate=result['full_material_gate']))
        print(json.dumps({k:result[k] for k in ['full_material_gate','simulation_admission','moment_errors']}))
    except Exception as error:
        failure=dict(schema='paper-flower-saved-independent-failure-v1',last_event=dict(stage),exception_type=type(error).__name__,exception_message=str(error),simulation_admission=False,native=identity)
        with (out/'saved-check-failure.json').open('x') as stream:
            stream.write(json.dumps(failure,indent=2)+'\n');stream.flush();os.fsync(stream.fileno())
        event(dict(stage='failed',exception_message=str(error)))
        raise
    finally:events.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ['folder','application','derived','assembly-reference','output']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--expected-native-sha',required=True)
    parser.add_argument('--preflight-only',action='store_true')
    args=parser.parse_args()
    run_check(args.folder,args.application,args.derived,args.assembly_reference,args.expected_native_sha,args.output,args.preflight_only)
