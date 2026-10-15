"""Independent bounded binary64 input gate, never an optimizer/derivative gate."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
from fractions import Fraction
import gmpy2

from rational_request_view import sha

TOLERANCE=gmpy2.mpq(1,10**6)


def text(q):return str(gmpy2.mpz(q.numerator))+'/'+str(gmpy2.mpz(q.denominator))


def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def vertices(packet):
    result=[]
    for row in packet['vertices']:
        if len(row)!=3:raise ValueError('INVALID_VERTEX')
        q=tuple(gmpy2.mpq(x) for x in row)
        if tuple(text(x) for x in q)!=tuple(row):raise ValueError('NONCANONICAL_COORDINATE')
        result.append(q)
    return result


def check_pair(authority,derived,length):
    """Recompute correspondence, tube bound and moments independently of claims."""
    from rounded_input import exact_moments
    length=gmpy2.mpq(length)
    if length<=0:raise ValueError('INVALID_INPUT_FIXED_LENGTH')
    if authority['faces']!=derived['faces']:raise ValueError('CONNECTIVITY_CHANGED')
    a,b=vertices(authority),vertices(derived)
    if len(a)!=len(b):raise ValueError('VERTEX_COUNT_CHANGED')
    r=gmpy2.mpq(0)
    for p,q in zip(a,b):
        for x,y in zip(p,q):
            if not math.isfinite(float(y)) or gmpy2.mpq(float(y))!=y:
                raise ValueError('DERIVED_COORDINATE_NOT_BINARY64')
            # mpq -> float goes through the current MPFR context and can
            # double-round near a binary64 midpoint. Fraction rounds directly
            # from exact integers, matching the frozen candidate conversion.
            nearest=float(Fraction(int(x.numerator),int(x.denominator)))
            if gmpy2.mpq(nearest)!=y:raise ValueError('GLOBAL_ROUNDING_MAPPING_CHANGED')
            r=max(r,abs(x-y))
    bound=gmpy2.mpq(0)
    for i,face in enumerate(authority['faces']):
        if len(face)!=3 or any(type(k)is not int or not 0<=k<len(a) for k in face):
            raise ValueError('INVALID_FACE: '+str(i))
        tri=[a[k] for k in face]; origin,u,v=tri
        normal=cross(tuple(u[k]-origin[k] for k in range(3)),tuple(v[k]-origin[k] for k in range(3)))
        # Independent sequential extrusion expansion, not candidate bound reuse.
        wx,wy,wz=[max(p[k] for p in tri)-min(p[k] for p in tri) for k in range(3)]
        ayz,axz,axy=[abs(x)/2 for x in normal]
        bound+=2*r*ayz+2*r*(axz+2*r*wz)+2*r*(axy+2*r*(wx+wy)+4*r*r)
    am=exact_moments(authority['vertices'],authority['faces'])
    bm=exact_moments(derived['vertices'],derived['faces'])
    delta=[y-x for x,y in zip(am,bm)]
    errors=[abs(d)/length**(3 if k==0 else 4 if k<4 else 5) for k,d in enumerate(delta)]
    within=(3*r*r<=TOLERANCE**2*length**2 and bound<=TOLERANCE*length**3 and all(e<=TOLERANCE for e in errors))
    return dict(maximum_linf_shift_exact=text(r),
        bidirectional_boundary_hausdorff_linf_upper_bound_exact=text(r),
        boundary_hausdorff_l2_upper_bound_squared_exact=text(3*r*r),
        conditional_symmetric_difference_volume_upper_bound_exact=text(bound),
        independent_authority_moments_exact=[text(x) for x in am],
        independent_derived_moments_exact=[text(x) for x in bm],
        independent_moment_delta_exact=[text(x) for x in delta],
        normalized_moment_errors_exact=[text(x) for x in errors],
        normalization_length_authority_exact=text(length),numerical_bounds_within_tolerance=bool(within),
        entity_bound_applied=False,solver_admission_granted=False)


def require_report_matches(actual,claimed):
    for key in ('maximum_linf_vertex_shift_exact','conditional_symmetric_difference_volume_upper_bound_exact',
                'independent_authority_moments_exact','independent_derived_moments_exact','independent_moment_delta_exact'):
        actual_key='maximum_linf_shift_exact' if key=='maximum_linf_vertex_shift_exact' else key
        a,b=actual[actual_key],claimed[key]
        if isinstance(a,list):
            match=isinstance(b,list) and len(a)==len(b) and all(gmpy2.mpq(x)==gmpy2.mpq(y) for x,y in zip(a,b))
        else:match=gmpy2.mpq(a)==gmpy2.mpq(b)
        if not match:raise ValueError('CANDIDATE_CLAIM_MISMATCH: '+key)


def save(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.flush();os.fsync(stream.fileno())
    fd=os.open(path.parent,os.O_RDONLY | os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)


def identity(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def resolve_pin_policy(policy,companion,expected_sources):
    if policy=='ORIGINAL_PIN_STRICT_V1':
        if companion is not None:raise ValueError('EXPLICIT_PIN_POLICY_REQUIRED')
        return None
    if policy!='BOUNDED_PIN_INTERSECTION_V1':raise ValueError('UNSUPPORTED_PIN_POLICY')
    if companion is None:raise ValueError('PIN_COMPANION_REQUIRED')
    from bounded_pin_attachment import load_companion
    return load_companion(companion,expected_source_hashes=expected_sources)


def main():
    parser=argparse.ArgumentParser()
    for name in ('admission','candidate','source-binding','pin','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--expected-native-sha256',required=True)
    parser.add_argument('--pin-policy',choices=('ORIGINAL_PIN_STRICT_V1','BOUNDED_PIN_INTERSECTION_V1'),default='ORIGINAL_PIN_STRICT_V1')
    parser.add_argument('--pin-companion',type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    began=time.monotonic();stage='SOURCE_IDENTITY'
    with (args.output/'events.jsonl').open('x') as stream:
        def event(name):
            nonlocal stage
            stage=name;stream.write(json.dumps(dict(stage=stage,elapsed_seconds=time.monotonic()-began))+'\n')
            stream.flush();os.fsync(stream.fileno());print(stage,flush=True)
        try:
            record=json.loads((args.admission/'admission.json').read_text())
            candidate=json.loads((args.candidate/'candidate.json').read_text())
            binding=json.loads(args.source_binding.read_text())
            admission_sha=sha(args.admission/'admission.json')
            if (record['status']!='INPUT_ADMITTED_FEASIBILITY_UNKNOWN' or
                binding.get('source_identity_verified')is not True or
                binding['admission_sha256']!=admission_sha or
                candidate['source_admission_sha256']!=admission_sha or
                candidate['source_problem_digest']!=record['problem_digest']):
                raise ValueError('AUTHORITY_BINDING_MISMATCH')
            if sha(args.pin)!=binding['frozen_source']['pin.obj']:raise ValueError('FIXED_PIN_SOURCE_MISMATCH')
            pinned={str(args.admission/'admission.json'):admission_sha,
                    str(args.candidate/'candidate.json'):sha(args.candidate/'candidate.json'),
                    str(args.source_binding):sha(args.source_binding),str(args.pin):sha(args.pin)}
            if args.pin_policy=='ORIGINAL_PIN_STRICT_V1' and args.pin_companion is not None:
                raise ValueError('EXPLICIT_PIN_POLICY_REQUIRED')
            if args.pin_policy=='BOUNDED_PIN_INTERSECTION_V1':
                if args.pin_companion is None:raise ValueError('PIN_COMPANION_REQUIRED')
                pinned[str(args.pin_companion)]=sha(args.pin_companion)
                companion_source=json.loads(args.pin_companion.read_text())
                for name,want in companion_source.get('files',{}).items():
                    path=args.pin_companion.parent/name
                    if path.parent.resolve()!=args.pin_companion.parent.resolve() or path.is_symlink() or sha(path)!=want:
                        raise ValueError('PIN_COMPANION_FILE_IDENTITY_MISMATCH')
                    pinned[str(path)]=want
            for name,want in candidate['files'].items():
                path=args.candidate/name
                if sha(path)!=want:raise ValueError('STALE_DERIVED_CANDIDATE: '+name)
                pinned[str(path)]=want
            pairs={}
            for name in ('outer','fixed','design'):
                source=args.admission/(name+'.rational.json')
                if sha(source)!=record['files'][source.name]:raise ValueError('STALE_AUTHORITY: '+name)
                pinned[str(source)]=record['files'][source.name]
                pairs[name]=(json.loads(source.read_text()),json.loads((args.candidate/source.name).read_text()))
                floating=json.loads((args.candidate/(name+'.binary64.json')).read_text())
                if floating['faces']!=pairs[name][1]['faces'] or [tuple(gmpy2.mpq(x) for x in row) for row in floating['vertices']]!=vertices(pairs[name][1]):
                    raise ValueError('REQUEST_FLOAT_VIEW_MISMATCH: '+name)
            av=vertices(pairs['outer'][0]);length=max(max(v[k] for v in av)-min(v[k] for v in av) for k in range(3))
            reverse={};reports={}
            for name,(source,derived) in pairs.items():
                event('INDEPENDENT_BOUNDS_AND_MOMENTS_'+name.upper())
                for p,q in zip(vertices(source),vertices(derived)):
                    if q in reverse and reverse[q]!=p:raise ValueError('GLOBAL_VERTEX_COLLISION')
                    reverse[q]=p
                report=check_pair(source,derived,length)
                claimed=candidate['regions'][name]
                require_report_matches(report,claimed)
                reports[name]=report;save(args.output/(name+'-independent-bounds.json'),report)
                if not report['numerical_bounds_within_tolerance']:raise ValueError('NUMERICAL_BOUND_OR_MOMENT_GATE_FAILED: '+name)
            import _make_it_stand as native
            if sha(native.__file__)!=args.expected_native_sha256:raise ValueError('NATIVE_IDENTITY_MISMATCH')
            from make_it_stand.verification.spin_problem import MeshSnapshot,SpinArtifactProblem
            snapshots={}
            for name in ('outer','fixed','design'):
                event('NATIVE_DERIVED_ENTITY_'+name.upper())
                value=pairs[name][1]
                snapshots[name]=MeshSnapshot.from_rationals(value['vertices'],value['faces'])
                save(args.output/(name+'-native-facts.json'),snapshots[name].exact_region().facts)
            event('EXACT_DERIVED_PARTITION')
            problem=SpinArtifactProblem(original=snapshots['outer'],
                fixed_regions=[('fixed-heart-and-pin',snapshots['fixed'],1.)],
                design_regions=[('distance-eroded-heart',snapshots['design'],1.)],
                length_unit='25mm',mass_unit='18.75g',
                target=dict(mode='SPINNING_TOP',axis=[0.,0.,1.],origin=[0.,0.,0.],support_points=[]),
                minimum_outer_wall=.015)
            save(args.output/'problem.json',problem.identity)
            event('EXACT_CONTAINMENT_AND_FIXED_PIN')
            import trimesh
            mesh=trimesh.load(args.pin,force='mesh',process=False)
            original_pin_packet=dict(vertices=[[text(gmpy2.mpq(float(x))) for x in row] for row in mesh.vertices],faces=mesh.faces.tolist())
            pin=MeshSnapshot.from_rationals(original_pin_packet['vertices'],original_pin_packet['faces'])
            expected_pin_sources=dict(pin=sha(args.pin),source_binding=sha(args.source_binding),
                original_admission=admission_sha,candidate=sha(args.candidate/'candidate.json'),native=sha(native.__file__),
                **{'authority-outer':record['files']['outer.rational.json'],
                   'authority-fixed':record['files']['fixed.rational.json'],
                   'candidate-fixed.rational.json':sha(args.candidate/'fixed.rational.json'),
                   'candidate-design.rational.json':sha(args.candidate/'design.rational.json')})
            pin_attachment=resolve_pin_policy(args.pin_policy,args.pin_companion,expected_pin_sources)
            extra_pin_relations={}
            if pin_attachment is not None:
                import numpy as np
                event('RECOMPUTED_BOUNDED_PIN_COMPANION_AND_FRESH_NATIVE_ROLE')
                if pin_attachment['original_pin_packet']!=original_pin_packet:
                    raise ValueError('PIN_COMPANION_ORIGINAL_LOADER_MISMATCH')
                value=pin_attachment['derived_packet']
                bounded_pin=MeshSnapshot.from_rationals(value['vertices'],value['faces'])
                authority_fixed=MeshSnapshot.from_rationals(pairs['fixed'][0]['vertices'],pairs['fixed'][0]['faces'])
                current_intersection=native.libigl_exact_boolean_construction_from_rationals(
                    original_pin_packet['vertices'],np.asarray(original_pin_packet['faces'],dtype=np.int32),
                    pairs['fixed'][1]['vertices'],np.asarray(pairs['fixed'][1]['faces'],dtype=np.int32),'intersection')
                if not current_intersection['succeeded']:raise ValueError('PIN_COMPANION_FRESH_INTERSECTION_FAILED')
                rebuilt=MeshSnapshot.from_rationals([[text(gmpy2.mpq(str(x))) for x in row] for row in current_intersection['vertices']],
                    np.asarray(current_intersection['faces'],dtype=np.int32).tolist())
                extra_pin_relations=dict(
                    original_pin_in_authority_fixed=native.relate_exact_material_regions(pin.exact_region(),authority_fixed.exact_region()),
                    original_pin_in_derived_fixed=native.relate_exact_material_regions(pin.exact_region(),snapshots['fixed'].exact_region()),
                    bounded_pin_in_original_pin=native.relate_exact_material_regions(bounded_pin.exact_region(),pin.exact_region()),
                    bounded_pin_equals_actual_intersection=native.relate_exact_material_regions(bounded_pin.exact_region(),rebuilt.exact_region()))
                save(args.output/'bounded-pin-fresh-relations.json',extra_pin_relations)
                if any(not r['completed'] for r in extra_pin_relations.values()) or \
                        not extra_pin_relations['original_pin_in_authority_fixed']['a_subset_b'] or \
                        not extra_pin_relations['bounded_pin_in_original_pin']['a_subset_b'] or \
                        not extra_pin_relations['bounded_pin_equals_actual_intersection']['equal']:
                    raise ValueError('PIN_COMPANION_FRESH_NATIVE_ROLE_FAILED')
                pin=bounded_pin
                save(args.output/'bounded-pin-consumption.json',dict(policy=args.pin_policy,
                    attachment_sha256=pin_attachment['attachment_sha256'],
                    recomputed_numerical_gate=pin_attachment['recomputed_numerical_gate'],
                    fresh_native_intersection_equality=True,original_pin_unchanged=True,
                    identity_semantics='EXPLICIT_BOUNDED_COMPANION_NOT_ORIGINAL_PIN_ENTITY'))
            relations=dict(design_in_outer=native.relate_exact_material_regions(snapshots['design'].exact_region(),snapshots['outer'].exact_region()),
                pin_in_fixed=native.relate_exact_material_regions(pin.exact_region(),snapshots['fixed'].exact_region()),
                pin_design=native.relate_exact_material_regions(pin.exact_region(),snapshots['design'].exact_region()))
            save(args.output/'relations.json',relations)
            if (not all(v['completed'] for v in relations.values()) or
                not relations['design_in_outer']['a_subset_b'] or
                not relations['pin_in_fixed']['a_subset_b'] or relations['pin_design']['positive_volume_overlap']):
                raise ValueError('DERIVED_CONTAINMENT_OR_PIN_ROLE_FAILED')
            event('EXACT_CLEARANCE')
            clearance=native.check_exact_boundary_clearance(snapshots['outer'].exact_region(),snapshots['design'].exact_region(),.015)
            save(args.output/'clearance.json',clearance)
            if not clearance['completed'] or clearance['at_least_minimum']is not True:
                raise ValueError('DERIVED_CLEARANCE_FAILED')
            for path,want in pinned.items():
                if sha(path)!=want:raise ValueError('SOURCE_CHANGED_DURING_CERTIFICATION: '+path)
            event('GEOMETRY_CHECKS_COMPLETE')
            # Previous authority native validity is source-bound, never inferred
            # from the tube formula; current derived validity was actually run.
            for report in reports.values():report['entity_bound_applied']=True
            evidence=dict(schema='bounded-derived-request-native-admission-v1',
                status='DERIVED_REQUEST_GEOMETRY_ADMITTED_FEASIBILITY_UNKNOWN',
                source_problem_digest=record['problem_digest'],derived_problem_digest=problem.input_digest,
                authority_admission_sha256=admission_sha,source_binding_sha256=sha(args.source_binding),
                candidate_sha256=sha(args.candidate/'candidate.json'),source_files=pinned,
                authority_validity='BOUND_TO_PREVIOUS_EXACT_NATIVE_ADMISSION',
                derived_validity='ACTUAL_NATIVE_CHECKS',bounds=reports,
                native_sha256=sha(native.__file__),checker_sha256=sha(__file__),
                geometry_admission_granted=True,solver_admission_granted=False,derivative_admission_granted=False,
                full_manufacturing_wall_thickness='UNKNOWN',optimization_invocations=0,simulation_invocations=0,
                normalization_length_authority_exact=text(length),tolerance_exact=text(TOLERANCE))
            if pin_attachment is not None:
                evidence.update(schema='bounded-derived-request-native-admission-v2',pin_policy=args.pin_policy,
                    bounded_pin_companion_sha256=pin_attachment['attachment_sha256'],
                    pin_role='EXPLICIT_BOUNDED_INTERSECTION_COMPANION',original_pin_source_unchanged=True)
            evidence['derived_request_identity']=identity(evidence)
            save(args.output/'admission.json',evidence)
            save(args.output/'terminal.json',dict(status=evidence['status'],completed=True,
                optimization_invocations=0,simulation_invocations=0,elapsed_seconds=time.monotonic()-began))
        except Exception as exc:
            save(args.output/'failure.json',dict(status='DERIVED_REQUEST_GEOMETRY_REJECTED',stage=stage,
                type=type(exc).__name__,message=str(exc),geometry_admission_granted=False,
                optimization_invocations=0,simulation_invocations=0,elapsed_seconds=time.monotonic()-began))
            raise


if __name__=='__main__':main()
