"""Frozen noodles preparation for the existing formal Cherchi solver.

The 1 mm seal requirement is implied by the measured frozen design gap. It
does not change the design domain; the original zero-wall request is retained.
"""
from dataclasses import asdict,replace
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
SEAL_WALL_M=.001


def compare_requests(old,new,normalized_wall):
    old=json.loads(json.dumps(old));new=json.loads(json.dumps(new))
    if old.keys()!=new.keys():raise ValueError('REQUEST_FIELDS_CHANGED')
    if old['minimum_wall_thickness']!=0. or new['minimum_wall_thickness']!=normalized_wall:
        raise ValueError('UNREGISTERED_WALL_CONTRACT')
    a=old.pop('fixed_shell_moments');b=new.pop('fixed_shell_moments')
    if len(a)!=10 or len(b)!=10 or not all(math.isfinite(x) for x in a+b):
        raise ValueError('INVALID_FIXED_MOMENTS')
    error=max(abs(x-y) for x,y in zip(a,b))
    if error>1e-10:raise ValueError('FIXED_MATERIAL_MOMENTS_CHANGED')
    old.pop('minimum_wall_thickness');new.pop('minimum_wall_thickness')
    if old!=new:raise ValueError('FROZEN_REQUEST_BEHAVIOR_CHANGED')
    return dict(equivalent=True,fixed_moment_roundoff=error,
        differences_allowed=['CERTIFIED_REDUNDANT_WALL_REQUIREMENT','FIXED_INTEGRATION_ROUNDOFF'])


def prepare(output,seed):
    import numpy as np
    import trimesh
    from make_it_stand import spin
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.faster_material import _prepare
    from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator,_json_value
    from finite_sample_summary import summarize
    if seed not in (0,1):raise ValueError('UNREGISTERED_CASE_SEED')
    load=lambda p:json.loads(p.read_text())
    request_path=HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    frozen=load(request_path)
    numeric=load(HERE/'noodles-formal-cherchi-derivatives-job-01/output/report.json')
    connections=load(HERE/'noodles-formal-cherchi-connection-job-01/output/report.json')
    if not summarize(numeric,connections)['finite_sample_evidence_passed']:
        raise ValueError('FINITE_PRODUCTION_EVIDENCE_NOT_PASSED')
    clearance=load(HERE/'noodles-saved-input-clearance-job-01/output/report.json')
    if (not clearance['measurement_completed'] or
        clearance['request_sha256']!=hashlib.sha256(request_path.read_bytes()).hexdigest() or
        Fraction(SEAL_WALL_M)**2>Fraction(clearance['squared_distance_exact'])):
        raise ValueError('SEAL_WALL_NOT_IMPLIED_BY_FROZEN_GEOMETRY')
    for name,digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:
            raise ValueError('FROZEN_INPUT_CHANGED')
    paths={Path(p).name:Path(p) for p in frozen['input_sources']}
    meshes={role:trimesh.load(paths[role+'-si-candidate.obj'],force='mesh',process=False)
            for role in ('outer','design','fixed')}
    spec_path=next(Path(p) for p in frozen['input_sources'] if Path(p).parent.name=='noodles-dense-v5')
    spec=load(spec_path)['spec']
    target=(np.asarray(spec['target_origin_model'])*spec['model_unit_m']).tolist()
    output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    problem=SpinArtifactProblem(original=(meshes['outer'].vertices,meshes['outer'].faces),
        fixed_regions=[('fixed-tubes-handle-tip',(meshes['fixed'].vertices,meshes['fixed'].faces),spec['density_kg_m3'])],
        design_regions=[('design',(meshes['design'].vertices,meshes['design'].faces),spec['density_kg_m3'])],
        length_unit='m',mass_unit='kg',minimum_outer_wall=SEAL_WALL_M,
        target=dict(mode='SPINNING_TOP',axis=spec['target_axis'],origin=target,support_points=[]))
    print('LOCKED_PROBLEM_READY',time.monotonic()-started,flush=True)
    old=frozen['requests'][seed]
    if old['settings']['deterministic_seed']!=seed:raise ValueError('SEED_ORDER_MISMATCH')
    evaluator=spin.SpinEvaluatorConfig(**dict(frozen['evaluator_settings'],algebraic_finite_map=True))
    solver=spin.SpinSolverConfig(**old['settings'])
    _,request,_=_prepare(problem,evaluator,solver)
    request=replace(request,center_of_mass_weight=old['center_of_mass_weight'],
                    spinning_top_objective=old['spinning_top_objective'])
    comparison=compare_requests(old,asdict(request),SEAL_WALL_M/problem.normalization['L'])
    (output/'prepared-request.json').write_text(json.dumps(asdict(request)))
    saved=load(HERE/'noodles-formal-cherchi-derivatives-job-01/output/point-000-material.json')
    mapping=InProcessCherchiEvaluator(request,evaluator,normalization=problem.normalization,
        design_vertices_physical=meshes['design'].vertices,
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',fixed_boundary_reuse=True)
    actual=mapping(saved['quadric_coefficients']);material=mapping.material(saved['quadric_coefficients'])
    (output/'same-call-material.json').write_text(json.dumps(material,default=_json_value))
    checks=dict(rational_vertices=material['variable_rational_vertices']==saved['variable_rational_vertices'],
        faces=material['variable_faces']==saved['variable_faces'],
        source_vertices=material['quadric_source_vertices']==saved['quadric_source_vertices'],
        source_faces=material['quadric_source_faces']==saved['quadric_source_faces'],
        jacobian_error=float(np.max(np.abs(np.asarray(actual['moment_jacobian'])-saved['moment_jacobian']))),
        moments_error=float(np.max(np.abs(np.asarray(actual['raw_moments'])-saved['raw_moments']))))
    passed=all(checks[k] for k in ('rational_vertices','faces','source_vertices','source_faces')) and checks['jacobian_error']<=1e-10 and checks['moments_error']<=1e-10
    report=dict(passed=passed,seed=seed,request_comparison=comparison,material_comparison=checks,
        input_digest=problem.input_digest,seal_wall_m=SEAL_WALL_M,original_request_unchanged=True,
        initialization='NO_OPTIMIZATION_DIAGNOSTIC_POINT_ONLY',optimization_invocations=0,
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report));print(json.dumps(report),flush=True)
    if not passed:raise ValueError('FORMAL_CASE_MAPPING_CHANGED')


if __name__=='__main__':prepare(Path(sys.argv[1]),int(sys.argv[2]))
