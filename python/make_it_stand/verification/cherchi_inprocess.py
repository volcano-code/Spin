"""Formal same-call in-process Cherchi geometry and moment evaluator.

The native Newton and existing sealed-material bridge consume this evaluator.
Geometry attachments are returned as fresh views of immutable JSON records.
A successful evaluation alone does not grant derivative or solve admission.
"""
from fractions import Fraction
from math import fsum
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import numpy as np
import time
import os
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


class InProcessCherchiEvaluator:
    def __init__(self,request,settings,*,branch_identity_policy='SOURCE_FEATURE_BRANCH_V2',stage_observer=None,
                 normalization=None,design_vertices_physical=None,
                 rational_arithmetic_backend='PYTHON_FRACTION_V1',fixed_boundary_reuse=False):
        import _cherchi_core_boolean as bridge
        if branch_identity_policy not in ('legacy','SOURCE_FEATURE_BRANCH_V2'):
            raise ValueError('UNSUPPORTED_BRANCH_IDENTITY_POLICY')
        self.branch_identity_policy=branch_identity_policy
        if rational_arithmetic_backend not in ('PYTHON_FRACTION_V1','NATIVE_GMP_SURFACE_V1'):
            raise ValueError('UNSUPPORTED_RATIONAL_ARITHMETIC')
        self.rational_arithmetic_backend=rational_arithmetic_backend
        self.arithmetic=None
        self.arithmetic_identity=None
        if rational_arithmetic_backend=='NATIVE_GMP_SURFACE_V1':
            try:
                import _cherchi_rational_moments as arithmetic
            except ImportError as error:
                raise RuntimeError('CHERCHI_NATIVE_RATIONAL_ARITHMETIC_UNAVAILABLE') from error
            algorithm='DIV_SURFACE_V1_REYNOLDS_V1_GMP_BINARY64_NEAREST_EVEN_V1'
            if getattr(arithmetic,'arithmetic_identity',None)!=algorithm:
                raise ValueError('UNSUPPORTED_NATIVE_RATIONAL_ARITHMETIC_IDENTITY')
            self.arithmetic=arithmetic
            self.arithmetic_identity=dict(backend=rational_arithmetic_backend,algorithm=algorithm,
                module_sha256=hashlib.sha256(Path(arithmetic.__file__).read_bytes()).hexdigest())
        self.request=request;self.settings=settings;self.bridge=bridge
        self.stage_observer=stage_observer
        self.call_counts={'source':0,'boolean':0,'integration':0,'derivative':0}
        self.stage_seconds={}
        self._materials={}
        self.native=spin._native_module()
        rotation=np.asarray(request.target_axes).T;origin=np.asarray(request.target_origin)
        if request.material_density<=0 or not np.allclose(rotation@rotation.T,np.eye(3),rtol=0,atol=1e-12) or abs(np.linalg.det(rotation)-1)>1e-12:
            raise ValueError('INVALID_DENSITY_OR_TARGET_FRAME')
        self.context_identity=digest(dict(request=asdict(request),settings=asdict(settings),
            native_sha256=hashlib.sha256(Path(self.native.__file__).read_bytes()).hexdigest(),
            bridge_sha256=hashlib.sha256(Path(bridge.__file__).read_bytes()).hexdigest(),
            adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            spin_sha256=hashlib.sha256(Path(spin.__file__).read_bytes()).hexdigest(),
            differential_sha256=hashlib.sha256(Path(__file__).with_name('cherchi_motion_differential.py').read_bytes()).hexdigest()))
        if branch_identity_policy!='legacy':
            self.context_identity=digest(dict(legacy_context=self.context_identity,branch_identity_policy=branch_identity_policy))
        if self.arithmetic_identity is not None:
            self.context_identity=digest(dict(base=self.context_identity,rational_arithmetic_identity=self.arithmetic_identity))
        self.backend_identity='cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:'+hashlib.sha256(Path(bridge.__file__).read_bytes()).hexdigest()
        self.target_design=(np.asarray(request.design_vertices)-origin)@rotation.T
        self.normalization=json.loads(json.dumps(normalization or dict(origin=(0.,0.,0.),L=1.)))
        self.length=float(self.normalization['L'])
        self.normalization_origin=np.asarray(self.normalization['origin'],dtype=float)
        if (not np.isfinite(self.length) or self.length<=0 or
            self.normalization_origin.shape!=(3,) or not np.isfinite(self.normalization_origin).all()):
            raise ValueError('INVALID_NORMALIZATION')
        self.design=np.array(request.design_vertices if design_vertices_physical is None
                             else design_vertices_physical,dtype=float,copy=True)
        if (self.design.shape!=np.asarray(request.design_vertices).shape or not np.isfinite(self.design).all()
            or not np.allclose((self.design-self.normalization_origin)/self.length,
                               request.design_vertices,atol=1e-14,rtol=0)):
            raise ValueError('PHYSICAL_DESIGN_IDENTITY_MISMATCH')
        self.design.setflags(write=False)
        self.target_design.setflags(write=False)
        self.normalization_origin.setflags(write=False)
        self.context_identity=digest(dict(base=self.context_identity,boolean_frame='PHYSICAL_INPUT_V1',
            design=self.design.tolist(),normalization=self.normalization))
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
        self._fixed_boundary_proof=None
        if type(fixed_boundary_reuse) is not bool:
            raise ValueError('fixed_boundary_reuse must be boolean')
        if fixed_boundary_reuse:
            if not hasattr(self.native,'verify_fixed_rational_boundary'):
                raise RuntimeError('VERIFIED_FIXED_BOUNDARY_UNAVAILABLE')
            _,reference=self._target_view(tuple(tuple(Fraction(float(x)) for x in row) for row in self.design))
            tokens=[[f'{x.numerator}/{x.denominator}' for x in row] for row in reference]
            self._fixed_boundary_proof=self._stage('fixed_boundary_proof',lambda:
                self.native.verify_fixed_rational_boundary(tokens,np.asarray(self.request.design_faces,dtype=np.int32)))
            self.context_identity=digest(dict(base=self.context_identity,
                fixed_boundary_proof_identity=self._fixed_boundary_proof.reference_identity))

    @property
    def fixed_boundary_proof(self):
        return self._fixed_boundary_proof

    def __call__(self,q):
        q=tuple(float(x) for x in q)
        if len(q)!=10 or not np.isfinite(q).all():raise ValueError('INVALID_QUADRIC_COEFFICIENTS')
        source=self._stage('source',lambda:self.native.evaluate_spin_it_faster_source_motion(np.asarray(q),self.target_design,self.triangulation))
        if not source['defined'] or not source['side_selection_defined']:
            raise ValueError('SOURCE_OR_SELECTOR_UNDEFINED:'+source['diagnostic'])
        axes=np.asarray(self.request.target_axes)
        sv=(np.asarray(source['vertices'])@axes.T+np.asarray(self.request.target_origin))*self.length+self.normalization_origin
        sf=np.asarray(source['faces'])
        operation='intersection' if source['interior_is_positive'] else 'difference'
        packet=self._stage('boolean',lambda:self.bridge.boolean(self.design.tolist(),[list(row) for row in self.request.design_faces],
            sv.tolist(),sf.tolist(),operation))
        self._validate_packet(packet)
        physical_vertices=np.asarray(packet['vertices']);faces=np.asarray(packet['faces'])
        physical_rational=tuple(tuple(exact_fraction(x) for x in row) for row in packet['exact_vertices'])
        faces,orientation=normalize_closed_component_orientation(physical_vertices,faces,rational_vertices=physical_rational)
        vertices,rational=self._target_view(physical_rational)
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
        physical_velocity=np.einsum('vkd,jk->vjd',velocity,axes)*self.length
        physical_moving=self._stage('intersection_motion',lambda:
            recover_construction_vertex_motion(records,self.design,sv,physical_velocity,physical_vertices))
        moving=np.einsum('vkd,kj->vjd',physical_moving/self.length,axes)
        prepared = None
        if self.arithmetic is not None and hasattr(self.arithmetic,'PreparedBoundary'):
            prepared = self._stage('prepare_arithmetic_boundary',lambda:
                self.arithmetic.PreparedBoundary(
                    [[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rational],
                    np.asarray(faces,dtype=np.int64)))
        jacobian=self._stage('derivative',lambda:
            self._differentiate(vertices,faces,moving,rational,prepared_boundary=prepared))
        variable=self._stage('integration',lambda:
            self._integrate(vertices,faces,rational,prepared_boundary=prepared))
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
        mapped=dict(raw_moments=np.asarray(variable)+self.fixed,variable_raw_moments=variable,moment_jacobian=jacobian,
            derivative_defined=True,derivative_diagnostic='INPROCESS_CHERCHI_PROVENANCE_DISCRETE_MAP_V1',
            geometry_hash=geometry,topology_hash=topology,evaluation_map_identity=self.context_identity,
            boolean_backend_identity=self.backend_identity,derivative_identity='cherchi-constructed-vertex-motion-v1')
        if self.arithmetic_identity is not None:
            mapped['rational_arithmetic_identity']=dict(self.arithmetic_identity)

        material=dict(mapped,quadric_coefficients=q,variable_vertices=vertices.tolist(),
            variable_faces=faces.tolist(),
            variable_rational_vertices=[[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rational],
            variable_vertices_model=(vertices@np.asarray(self.request.target_axes).T+np.asarray(self.request.target_origin)).tolist(),
            quadric_source_vertices=sv.tolist(),quadric_source_faces=sf.tolist(),
            quadric_source_identity=digest(dict(vertices=sv.tolist(),faces=sf.tolist(),coefficients=q)),
            boolean_operation='intersection' if operation=='intersection' else 'subtraction',
            material_proof=dict(schema='cherchi-same-call-material-v1',
                exact_vertices=packet['exact_vertices'],construction_records=packet['construction_records'],
                resolved_construction_records=records,orientation=orientation))
        self._materials[q]=self._stage('material_attachment_serialization',lambda:
            json.dumps(material,default=_json_value,sort_keys=True,allow_nan=False))
        return mapped

    def _integrate(self,vertices,faces,rational,*,prepared_boundary=None):
        if self.arithmetic is None:
            return integrate_rational_polyhedral_raw_moments(vertices,faces,density=self.request.material_density,
                rational_vertices=rational)
        if prepared_boundary is not None:
            return prepared_boundary.integrate_moments(self.request.material_density)
        tokens=[[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rational]
        return self.arithmetic.integrate_moments(tokens,np.asarray(faces,dtype=np.int64),self.request.material_density)

    def _differentiate(self,vertices,faces,moving,rational,*,prepared_boundary=None):
        if self.arithmetic is None:
            return differentiate_polyhedral_raw_moments(vertices,faces,moving,density=self.request.material_density,
                rational_vertices=rational)
        tokens=([[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rational]
                if prepared_boundary is None else None)
        indexed_faces=np.asarray(faces,dtype=np.int64)
        face_motion=moving[indexed_faces]
        if hasattr(self.arithmetic, 'transport_coefficients_selected'):
            # No tolerance: arbitrarily small nonzero motion remains active.
            # The native entry validates the entire original closed boundary.
            selected=np.flatnonzero(np.any(face_motion!=0.,axis=(1,2,3)))
            coefficients=(prepared_boundary.transport_coefficients_selected(
                self.request.material_density,selected.tolist()) if prepared_boundary is not None else
                self.arithmetic.transport_coefficients_selected(
                    tokens,indexed_faces,self.request.material_density,selected.tolist()))
            terms=np.einsum('fjmc,fjcd->fmd',coefficients,face_motion[selected])
        else:
            # Compatibility with frozen earlier modules, using the same GMP algorithm.
            coefficients=self.arithmetic.transport_coefficients(tokens,indexed_faces,self.request.material_density)
            terms=np.einsum('fjmc,fjcd->fmd',coefficients,face_motion)
        jacobian=np.asarray([[fsum(terms[:,m,d]) for d in range(moving.shape[2])] for m in range(10)])
        if not np.isfinite(jacobian).all():raise ValueError('NATIVE_RATIONAL_DIFFERENTIAL_NONFINITE')
        return jacobian

    def _stage(self,name,operation):
        start=time.monotonic()
        self.call_counts[name]=self.call_counts.get(name,0)+1
        try:
            return operation() if self.stage_observer is None else self.stage_observer(name,operation)
        finally:
            self.stage_seconds[name]=self.stage_seconds.get(name,0.)+time.monotonic()-start

    def _target_view(self,physical_rational):
        """Exact affine change of coordinates, preserving the input entity."""
        origin=[Fraction(float(x)) for x in self.normalization_origin]
        target_origin=[Fraction(float(x)) for x in self.request.target_origin]
        axes=[[Fraction(float(x)) for x in row] for row in self.request.target_axes]
        length=Fraction(self.length)
        rows=[]
        for physical in physical_rational:
            model=[(physical[i]-origin[i])/length-target_origin[i] for i in range(3)]
            rows.append(tuple(sum(model[i]*axes[i][j] for i in range(3)) for j in range(3)))
        return np.asarray([[float(x) for x in row] for row in rows],dtype=float).reshape(-1,3),tuple(rows)

    @staticmethod
    def _validate_packet(packet):
        if packet.get('backend')!='cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0-inprocess-v1':
            raise ValueError('CHERCHI_BACKEND_IDENTITY_MISMATCH')

    def material(self,coefficients):
        key=tuple(float(x) for x in coefficients)
        if key not in self._materials:raise ValueError('MATERIAL_NOT_EVALUATED')
        result=json.loads(self._materials[key])
        result['quadric_coefficients']=tuple(result['quadric_coefficients'])
        return result

    def validate_same_call_material(self, material, coefficients, request):
        """Bind trial moments and derivative to this producer's immutable packet."""
        if request != self.request:
            raise ValueError('SAME_CALL_REQUEST_MISMATCH')
        key = tuple(float(x) for x in coefficients)
        if key not in self._materials:
            raise ValueError('SAME_CALL_MATERIAL_NOT_EVALUATED')
        encoded = json.dumps(material, default=_json_value, sort_keys=True, allow_nan=False)
        if encoded != self._materials[key]:
            raise ValueError('SAME_CALL_MATERIAL_MISMATCH')

    def geometry(self,coefficients,context):
        """Callback for the existing strict solve-and-seal interface."""
        design=np.asarray(context['design_vertices'],dtype=np.float64)
        if (design.shape!=self.target_design.shape or design.tobytes()!=self.target_design.tobytes()
            or not np.array_equal(context['design_faces'],self.request.design_faces)
            or context['material_density']!=self.request.material_density):
            raise ValueError('CHERCHI_LOCKED_CONTEXT_MISMATCH')
        self(coefficients)
        return self.material(coefficients)

    def complement(self,retained,context=None,*,packet_sink=None):
        """Opposite-side Boolean using the accepted source, with no source reconstruction."""
        saved=self.material(retained['quadric_coefficients'])
        fields=('geometry_hash','quadric_source_identity','quadric_source_vertices','quadric_source_faces')
        if any(retained.get(field)!=saved[field] for field in fields):
            raise ValueError('TERMINAL_SOURCE_IDENTITY_MISMATCH')
        return self.complement_of_saved_source(saved,packet_sink=packet_sink)

    def complement_of_saved_source(self, saved, *, packet_sink=None):
        """Postprocessing only; this does not grant a same-call solver evaluation.

        The recovery caller must separately bind a verified accepted checkpoint
        to this request. Source vertices are used verbatim, never reconstructed
        from coefficients or supplied as Newton initialization.
        """
        if packet_sink is not None and not callable(packet_sink):
            raise ValueError('INVALID_COMPLEMENT_PACKET_SINK')
        saved=json.loads(json.dumps(saved,default=_json_value,allow_nan=False))
        saved['quadric_coefficients']=tuple(saved['quadric_coefficients'])
        if saved.get('quadric_source_identity') != digest(dict(
                vertices=saved['quadric_source_vertices'], faces=saved['quadric_source_faces'],
                coefficients=saved['quadric_coefficients'])):
            raise ValueError('TERMINAL_SOURCE_IDENTITY_MISMATCH')
        if saved.get('boolean_operation') not in ('intersection', 'subtraction'):
            raise ValueError('INVALID_RETAINED_BOOLEAN_OPERATION')
        operation='difference' if saved['boolean_operation']=='intersection' else 'intersection'
        sv=saved['quadric_source_vertices'];sf=saved['quadric_source_faces']
        packet=self._stage('boolean',lambda:self.bridge.boolean(self.design.tolist(),
            [list(row) for row in self.request.design_faces],sv,sf,operation))
        self._validate_packet(packet)
        if packet_sink is not None:
            packet_sink(packet)
        physical_vertices=np.asarray(packet['vertices'],dtype=float);faces=np.asarray(packet['faces'],dtype=int)
        physical_rational=tuple(tuple(exact_fraction(x) for x in row) for row in packet['exact_vertices'])
        vertices,rational=self._target_view(physical_rational)
        if len(vertices):
            faces,orientation=normalize_closed_component_orientation(physical_vertices,faces,rational_vertices=physical_rational)
            variable=self._stage('integration',lambda:self._integrate(vertices,faces,rational))
        else:
            vertices=np.empty((0,3));faces=np.empty((0,3),dtype=int);variable=np.zeros(10)
        geometry=digest(dict(vertices=packet['exact_vertices'],faces=faces.tolist()))
        return dict(variable_vertices=vertices.tolist(),variable_faces=faces.tolist(),
            variable_rational_vertices=[[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rational],
            variable_vertices_model=(vertices@np.asarray(self.request.target_axes).T+np.asarray(self.request.target_origin)).tolist(),
            variable_raw_moments=variable,moment_jacobian=np.zeros((10,10)),
            derivative_defined=False,derivative_identity='COMPLEMENT_NO_DERIVATIVE',
            geometry_hash=geometry,topology_hash=digest(dict(operation=operation,faces=faces.tolist())),
            evaluation_map_identity=self.context_identity,boolean_backend_identity=self.backend_identity,
            quadric_coefficients=saved['quadric_coefficients'],quadric_source_vertices=sv,quadric_source_faces=sf,
            quadric_source_identity=saved['quadric_source_identity'],
            boolean_operation='intersection' if operation=='intersection' else 'subtraction',
            material_proof=dict(schema='cherchi-same-call-material-v1',exact_vertices=packet['exact_vertices'],
                construction_records=packet['construction_records']))


def _json_value(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(f'Unsupported material record value: {type(value).__name__}')


def solve_cherchi_inprocess_and_seal_faster(problem, *, evaluator_settings,
        solver_settings, output_directory, center_of_mass_weight=0.0,
        spinning_top_objective='INERTIA_AND_LOWERING', event_sink=None,
        rational_arithmetic_backend='PYTHON_FRACTION_V1',solver_boundary_profile='EXACT_NEF_REGION_V1',
        durability_profile='PER_EVENT_FSYNC_V1',
        validation_profile='INDEPENDENT_EACH_EVALUATION_V1', input_validation_profile='RECHECK_INPUT_V1',
        terminal_region_profile='INDEXED_BOUNDARY_REBUILD_V1'):
    """Seeded native Newton with production Cherchi, durable records and sealing.

    Uses the existing locked-problem and exact material bridge. No explicit
    initial quadrics are accepted: starts come from the registered request seed.
    The directory is exclusive; a terminated process leaves an incomplete prefix.
    """
    from .faster_material import (_prepare_faster_problem, solve_cherchi_geometry_and_seal_faster,
                                 _validate_terminal_region_profile)
    from .cherchi_run_journal import RunJournal
    if durability_profile not in ('PER_EVENT_FSYNC_V1','EVALUATION_CHECKPOINT_V1'):
        raise ValueError('unknown run durability profile')
    if solver_boundary_profile not in ('EXACT_NEF_REGION_V1','EXACT_RATIONAL_BOUNDARY_V1','EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'):
        raise ValueError('unknown solver boundary profile')
    if validation_profile not in ('INDEPENDENT_EACH_EVALUATION_V1', 'LAYERED_SAME_CALL_V1'):
        raise ValueError('unknown material validation profile')
    _validate_terminal_region_profile(terminal_region_profile)
    prepared_problem = _prepare_faster_problem(problem, evaluator_settings, solver_settings,
        input_validation_profile=input_validation_profile)
    problem, request = prepared_problem.problem, prepared_problem.request
    request=replace(request,center_of_mass_weight=center_of_mass_weight,
                    spinning_top_objective=spinning_top_objective)
    output=Path(output_directory)
    output.mkdir(parents=True,exist_ok=False)
    with RunJournal(output,durability_profile=durability_profile,json_value=_json_value) as journal:
        publish=journal.publish
        emit=journal.emit
        def stage(name,operation):
            start=time.monotonic_ns()
            emit(dict(kind='STAGE_STARTED',stage=name,monotonic_ns=start))
            try:
                value=operation()
            except Exception as exc:
                emit(dict(kind='STAGE_FAILED',stage=name,error=f'{type(exc).__name__}: {exc}',
                          monotonic_ns=time.monotonic_ns(),elapsed_ns=time.monotonic_ns()-start))
                raise
            emit(dict(kind='STAGE_RETURNED',stage=name,monotonic_ns=time.monotonic_ns(),
                      elapsed_ns=time.monotonic_ns()-start))
            return value
        evaluator=InProcessCherchiEvaluator(request,evaluator_settings,stage_observer=stage,
            normalization=problem.normalization,design_vertices_physical=problem.design_regions[0].mesh.arrays()[0],
            rational_arithmetic_backend=rational_arithmetic_backend,
            fixed_boundary_reuse=solver_boundary_profile in ('EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'))
        if solver_boundary_profile!='EXACT_NEF_REGION_V1':
            evaluator.context_identity=digest(dict(base=evaluator.context_identity,
                solver_boundary_profile=solver_boundary_profile,
                verifier_sha256=hashlib.sha256(Path(__file__).with_name('faster_material.py').read_bytes()).hexdigest(),
                snapshot_sha256=hashlib.sha256(Path(__file__).with_name('spin_problem.py').read_bytes()).hexdigest()))
        if validation_profile != 'INDEPENDENT_EACH_EVALUATION_V1':
            evaluator.context_identity = digest(dict(base=evaluator.context_identity,
                validation_profile=validation_profile,
                verifier_sha256=hashlib.sha256(Path(__file__).with_name('faster_material.py').read_bytes()).hexdigest()))
        configuration=dict(schema='cherchi-inprocess-formal-run-v1',request=asdict(request),
            evaluator=asdict(evaluator_settings),evaluation_map_identity=evaluator.context_identity,
            boolean_backend_identity=evaluator.backend_identity,initialization='REQUEST_SEED_ONLY',complete=False)
        if evaluator.arithmetic_identity is not None:
            configuration.update(rational_arithmetic_backend=rational_arithmetic_backend,
                rational_arithmetic_identity=dict(evaluator.arithmetic_identity))
        if solver_boundary_profile!='EXACT_NEF_REGION_V1':
            configuration['solver_boundary_profile']=solver_boundary_profile
        if evaluator.fixed_boundary_proof is not None:
            configuration['fixed_boundary_proof_identity']=evaluator.fixed_boundary_proof.reference_identity
        if durability_profile!='PER_EVENT_FSYNC_V1':
            configuration['durability_profile']=durability_profile
        if validation_profile != 'INDEPENDENT_EACH_EVALUATION_V1':
            configuration['validation_profile'] = validation_profile
        if input_validation_profile != 'RECHECK_INPUT_V1':
            configuration['input_validation_profile'] = input_validation_profile
            configuration['input_validation_problem_digest'] = problem.input_digest
        if terminal_region_profile != 'INDEXED_BOUNDARY_REBUILD_V1':
            configuration['terminal_region_profile'] = terminal_region_profile
        publish('configuration.json',configuration)
        count=0
        references={}
        def geometry(q,context):
            nonlocal count
            count+=1
            key=tuple(float(x) for x in q)
            try:
                material=stage('production_geometry_evaluation',lambda:evaluator.geometry(key,context))
            except Exception as exc:
                emit(dict(kind='MAP_REJECTED',evaluation=count,coefficients=key,
                          error=f'{type(exc).__name__}: {exc}'))
                raise
            name=f'evaluation-{count:06d}.json'
            def checkpoint():
                publish(name,material)
                references[key]=dict(evaluation_file=name,geometry_hash=material['geometry_hash'],
                                     sha256=hashlib.sha256((output/name).read_bytes()).hexdigest())
                emit(dict(kind='MAP_SAVED',evaluation=count,coefficients=key,**references[key]))
            stage('evaluation_checkpoint',checkpoint)
            return material
        def sink(row):
            def record():
                emit(dict(kind='NEWTON_EVENT',event=row))
                if event_sink is not None:event_sink(row)
            stage('newton_event_record',record)
        def terminal_checkpoint(phase, payload):
            if phase == 'NATIVE_RETURNED':
                publish('native-result.json', dict(schema='cherchi-numerical-terminal-v1',
                    native_run=payload, configuration_sha256=hashlib.sha256(
                        (output/'configuration.json').read_bytes()).hexdigest(),
                    material_audit_complete=False))
                if payload.get('status') == 'CONVERGED':
                    q = tuple(float(x) for x in payload['quadric_coefficients'])
                    if q not in references:
                        raise ValueError('TERMINAL_MATERIAL_NOT_IN_DURABLE_LOG')
                    publish('terminal-material-reference.json', dict(coefficients=q, **references[q]))
            elif phase == 'COMPLEMENT_RETURNED':
                publish('terminal-complement.json', payload)
                publish('terminal-complement-reference.json', dict(
                    schema='cherchi-terminal-complement-reference-v1',
                    evaluation_file='terminal-complement.json',
                    sha256=hashlib.sha256((output/'terminal-complement.json').read_bytes()).hexdigest(),
                    retained_reference_sha256=hashlib.sha256(
                        (output/'terminal-material-reference.json').read_bytes()).hexdigest(),
                    configuration_sha256=hashlib.sha256((output/'configuration.json').read_bytes()).hexdigest(),
                    material_audit_complete=False))
            else:
                raise ValueError('UNKNOWN_TERMINAL_CHECKPOINT_PHASE')
            emit(dict(kind='TERMINAL_CHECKPOINT_SAVED', phase=phase))
        def complement(retained, context):
            return evaluator.complement(retained, context,
                packet_sink=lambda packet: publish('terminal-complement-raw.json', packet))
        result=solve_cherchi_geometry_and_seal_faster(problem,geometry_evaluator=geometry,
            complement_evaluator=complement,initial_quadrics=(),
            evaluator_settings=evaluator_settings,solver_settings=solver_settings,
            center_of_mass_weight=center_of_mass_weight,spinning_top_objective=spinning_top_objective,
            event_sink=sink,stage_observer=stage,require_same_call_terminal=True,
            solver_boundary_profile=solver_boundary_profile,
            validation_profile=validation_profile,
            material_validator=evaluator.validate_same_call_material,
            terminal_checkpoint_sink=terminal_checkpoint, _prepared_problem=prepared_problem,
            input_validation_profile=input_validation_profile, terminal_region_profile=terminal_region_profile)
        run=result.identity['native_run']
        from .faster_timing import summarize_solver_timing
        events=[json.loads(line) for line in (output/'events.jsonl').read_text().splitlines()]
        publish('timing-scope.json',summarize_solver_timing(events))
        publish('result.json',dict(status=result.status,identity=result.identity,
            scientific_digest=result.scientific_digest,call_counts=evaluator.call_counts,
            stage_seconds=evaluator.stage_seconds,durability_profile=durability_profile))
        emit(dict(kind='RUN_TERMINATED',status=result.status,evaluation_count=count))
        publish('complete.json',dict(status=result.status,complete=True,evaluation_count=count))
        if result.status=='CONVERGED':
            from .cherchi_run_evidence import verify_cherchi_run_evidence
            publish('record-verification.json',verify_cherchi_run_evidence(output))
    # Timing summary is outside all scientific event/material identities.
    with (output/'recording-costs.json').open('x') as stream:
        json.dump(dict(schema='cherchi-run-recording-costs-v1',durability_profile=durability_profile,
            metrics=journal.metrics,scope='Record I/O through journal close; summary write excluded'),stream,
            sort_keys=True,allow_nan=False)
        stream.flush();os.fsync(stream.fileno())
    descriptor=os.open(output,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(descriptor)
    finally:os.close(descriptor)
    return result
