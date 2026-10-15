"""Same-call in-process Cherchi moment map for numerical diagnostics.

This adapter has no sealing authority. Production admission still requires
the registered geometry, derivative and complete-solve checks.
"""
from fractions import Fraction
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_motion_differential import (
    normalize_closed_component_orientation,recover_construction_vertex_motion,
    differentiate_polyhedral_raw_moments,integrate_rational_polyhedral_raw_moments)


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def exact_fraction(token):
    numerator,denominator=token.split('/')
    def integer(text):return int(text,2 if text.lstrip('-').startswith('0b') else 10)
    return Fraction(integer(numerator),integer(denominator))


def construction_constraint_key(record):
    """Typed defining constraints, independent of equivalent operand ordering.

    LPI is an unoriented line intersecting an unoriented plane; TPI is the
    intersection of three unoriented planes. This key does not reorder any
    actual derivative operands or assert that the constraints are nonsingular.
    """
    operands=[]
    for operand in record['operands']:
        candidates=operand['source_vertices']
        if operand['status']!='UNIQUE' or len(candidates)!=1:
            raise ValueError('CONSTRUCTION_CONSTRAINT_CORRESPONDENCE_UNAVAILABLE')
        candidate=candidates[0];operands.append((candidate['input_object'],candidate['input_vertex']))
    kind=record['point_type']
    if kind=='EXPLICIT3D' and len(operands)==1:return (kind,tuple(operands))
    if kind=='LPI' and len(operands)==5:
        return (kind,tuple(sorted(operands[:2])),tuple(sorted(operands[2:])))
    if kind=='TPI' and len(operands)==9:
        return (kind,tuple(sorted(tuple(sorted(operands[i:i+3])) for i in (0,3,6))))
    raise ValueError('UNSUPPORTED_CONSTRUCTION_CONSTRAINT_SCHEMA')


def replay_lpi_constraint(record,operands):
    """Independent exact coordinate replay, not motion or solver admission."""
    construction_constraint_key(record)
    if record['point_type']!='LPI':raise ValueError('ONLY_LPI_CONSTRAINT_REPLAY_SUPPORTED')
    points=[]
    for operand in record['operands']:
        candidate=operand['source_vertices'][0]
        points.append(tuple(Fraction(x) for x in operands[candidate['input_object']][candidate['input_vertex']]))
    from exact_planar_patches import cross,subtract,dot
    a,b,c,d,e=points;normal=cross(subtract(d,c),subtract(e,c));direction=subtract(b,a)
    denominator=dot(normal,direction)
    if denominator==0:raise ValueError('SINGULAR_LPI_CONSTRAINT')
    t=dot(normal,subtract(c,a))/denominator
    return tuple(x+t*delta for x,delta in zip(a,direction))


def feature_branch_identity(operation,faces,source_faces,surface_type,cap_count,source_edges,records):
    """Discrete feature fingerprint, NOT a neighborhood/embedding certificate.

    Operand order is preserved; only output storage indices and cyclic face
    start indices are canonicalized. Interpolation parameters are geometry.
    Ambiguous or duplicate features cannot grant a correspondence identity.
    """
    keys=[]
    for record in records:
        operands=[]
        for operand in record['operands']:
            candidates=operand['source_vertices']
            if operand['status']!='UNIQUE' or len(candidates)!=1:
                raise ValueError('BRANCH_FEATURE_CORRESPONDENCE_UNAVAILABLE')
            candidates=candidates[0]
            operands.append((candidates['input_object'],candidates['input_vertex']))
        keys.append(json.dumps([record['point_type'],operands],separators=(',',':')))
    if len(set(keys))!=len(keys):raise ValueError('DUPLICATE_BRANCH_VERTEX_FEATURE')
    triangles=[]
    for face in faces:
        if len(face)!=3 or len(set(map(int,face)))!=3 or any(int(i)<0 or int(i)>=len(keys) for i in face):
            raise ValueError('INVALID_BRANCH_FACE')
        tri=tuple(keys[int(i)] for i in face)
        triangles.append(min(tri,tri[1:]+tri[:1],tri[2:]+tri[:2]))
    return digest(dict(policy='SOURCE_FEATURE_BRANCH_V2',operation=operation,
        oriented_material_triangles=sorted(triangles),source_faces=source_faces,
        surface_type=surface_type,cap_count=cap_count,source_edges=source_edges))


class CherchiCoreMap:
    def __init__(self,request,settings,*,branch_identity_policy='legacy'):
        import _cherchi_core_boolean as bridge
        if branch_identity_policy not in ('legacy','SOURCE_FEATURE_BRANCH_V2'):
            raise ValueError('UNSUPPORTED_BRANCH_IDENTITY_POLICY')
        self.branch_identity_policy=branch_identity_policy
        self.request=request;self.settings=settings;self.bridge=bridge
        self.native=spin._native_module()
        rotation=np.asarray(request.target_axes).T;origin=np.asarray(request.target_origin)
        if request.material_density<=0 or not np.allclose(rotation@rotation.T,np.eye(3),rtol=0,atol=1e-12) or abs(np.linalg.det(rotation)-1)>1e-12:
            raise ValueError('INVALID_DENSITY_OR_TARGET_FRAME')
        self.context_identity=digest(dict(request=asdict(request),settings=asdict(settings),
            native_sha256=hashlib.sha256(Path(self.native.__file__).read_bytes()).hexdigest(),
            bridge_sha256=hashlib.sha256(Path(bridge.__file__).read_bytes()).hexdigest(),
            adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
        if branch_identity_policy!='legacy':
            self.context_identity=digest(dict(legacy_context=self.context_identity,branch_identity_policy=branch_identity_policy))
        self.backend_identity='cherchi-inprocess-'+hashlib.sha256(Path(bridge.__file__).read_bytes()).hexdigest()
        self.design=(np.asarray(request.design_vertices)-origin)@rotation.T
        r=np.asarray(request.fixed_shell_moments);mass=r[0];first=r[1:4]
        second=np.array(((r[4],r[9],r[8]),(r[9],r[5],r[7]),(r[8],r[7],r[6])))
        f=rotation@(first-mass*origin)
        s=rotation@(second-np.outer(origin,first)-np.outer(first,origin)+mass*np.outer(origin,origin))@rotation.T
        self.fixed=np.asarray((mass,*f,s[0,0],s[1,1],s[2,2],s[1,2],s[0,2],s[0,1]))
        self.triangulation=self.native.SpinTriangulationSettings()
        for name in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                     'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
            setattr(self.triangulation,name,getattr(settings,name))
        self.last=None

    def __call__(self,q):
        source=self.native.evaluate_spin_it_faster_source_motion(np.asarray(q),self.design,self.triangulation)
        if not source['defined'] or not source['side_selection_defined']:
            raise ValueError('SOURCE_OR_SELECTOR_UNDEFINED:'+source['diagnostic'])
        sv=np.asarray(source['vertices']);sf=np.asarray(source['faces'])
        operation='intersection' if source['interior_is_positive'] else 'difference'
        packet=self.bridge.boolean(self.design.tolist(),[list(row) for row in self.request.design_faces],
            sv.tolist(),sf.tolist(),operation)
        vertices=np.asarray(packet['vertices']);faces=np.asarray(packet['faces'])
        rational=tuple(tuple(exact_fraction(x) for x in row) for row in packet['exact_vertices'])
        faces,orientation=normalize_closed_component_orientation(vertices,faces,rational_vertices=rational)
        lookup={}
        for object_index,points in enumerate((self.design,sv)):
            for vertex_index,point in enumerate(points):
                key=tuple(Fraction(float(x)) for x in point)
                lookup.setdefault(key,[]).append(dict(input_object=object_index,input_vertex=vertex_index))
        records=[]
        for item in packet['construction_records']:
            operands=[];positions=[]
            for tokens in item['operand_coordinates_exact']:
                exact=tuple(exact_fraction(x) for x in tokens)
                candidates=lookup.get(exact,[])
                operands.append(dict(source_vertices=candidates,status='UNIQUE' if len(candidates)==1 else
                    'UNMATCHED' if not candidates else 'AMBIGUOUS'))
                positions.append([float(x) for x in exact])
            records.append(dict(vertex_index=item['vertex_index'],point_type=item['point_type'],
                parameter_t=item['parameter_t'],operands=operands,operand_coordinates=positions))
        velocity=np.asarray(source['vertex_velocity']).reshape(len(sv),3,10)
        moving=recover_construction_vertex_motion(records,self.design,sv,velocity,vertices)
        jacobian=differentiate_polyhedral_raw_moments(vertices,faces,moving,density=self.request.material_density,
            rational_vertices=rational)
        variable=integrate_rational_polyhedral_raw_moments(vertices,faces,density=self.request.material_density,
            rational_vertices=rational)
        geometry=digest(dict(vertices=packet['exact_vertices'],faces=faces.tolist()))
        topology=digest(dict(operation=operation,faces=faces.tolist(),source_faces=sf.tolist(),
            surface_type=source['surface_type'],cap_count=source['cap_face_count'],
            source_edges=np.asarray(source['source_normal_edges']).tolist(),
            records=[dict(point_type=r['point_type'],operands=r['operands'],parameter_t=r['parameter_t']) for r in records]))
        if self.branch_identity_policy=='SOURCE_FEATURE_BRANCH_V2':
            topology=feature_branch_identity(operation,faces.tolist(),sf.tolist(),
                source['surface_type'],source['cap_face_count'],
                np.asarray(source['source_normal_edges']).tolist(),records)
        self.last=dict(vertices=vertices,faces=faces,exact_vertices=rational,source=source,orientation=orientation,
            operation=operation,geometry_hash=geometry,topology_hash=topology,
            provenance_schema='same-call-cherchi-construction-v1',
            original_construction_records=packet['construction_records'],
            resolved_construction_records=records)
        return dict(raw_moments=np.asarray(variable)+self.fixed,variable_raw_moments=variable,moment_jacobian=jacobian,
            derivative_defined=True,derivative_diagnostic='INPROCESS_CHERCHI_PROVENANCE_DISCRETE_MAP_DIAGNOSTIC',
            geometry_hash=geometry,topology_hash=topology,evaluation_map_identity=self.context_identity,
            boolean_backend_identity=self.backend_identity,derivative_identity='cherchi-constructed-vertex-motion-v1')
