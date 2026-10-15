"""Cherchi's own 81-point production check on the frozen noodles input.

Reuses registered coefficients/scales, not CGAL's pass flag. Saves source
construction metadata and same-call material for later finite classification.
"""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
import trimesh
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from make_it_stand.verification.spin_problem import rational_solver_candidate_snapshot
from independent_tetra_moments import integrate
from probe_formal_cherchi_derivatives import weights

HERE = Path(__file__).resolve().parent


def main(output):
    frozen_path = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    frozen = json.loads(frozen_path.read_text())
    for name, digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError('FROZEN_INPUT_CHANGED')
    old_manifest_path = HERE/'noodles-algebraic-production-moments-job-01/output/manifest.json'
    old = json.loads(old_manifest_path.read_text())
    if old['request_sha256'] != hashlib.sha256(frozen_path.read_bytes()).hexdigest():
        raise ValueError('REGISTERED_SAMPLE_REQUEST_MISMATCH')
    settings = dict(frozen['evaluator_settings'], algebraic_finite_map=True)
    if settings != old['settings'] or len(old['samples']) != 81 or len(old['layouts']) != 30:
        raise ValueError('REGISTERED_MAPPING_OR_LAYOUT_MISMATCH')
    samples = np.asarray(old['samples'])
    if samples.shape != (81,10) or len({tuple(q) for q in samples}) != 81:
        raise ValueError('INVALID_REGISTERED_SAMPLES')
    for item in old['layouts']:
        direction = item['direction']; indices = item['indices']
        if item['scale'] not in (3e-6,1.5e-6,.75e-6) or len(indices) != 5:
            raise ValueError('INVALID_REGISTERED_SCALE')
        for index, position in zip(indices,item['positions']):
            diff = samples[index]-samples[0]
            if float(diff[direction]) != position or np.any(np.delete(diff,direction) != 0):
                raise ValueError('LAYOUT_SAMPLE_MISMATCH')
    payload = dict(frozen['requests'][0])
    payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
    request = spin.SpinItFasterRequest(**payload)
    physical_path = next(Path(p) for p in frozen['input_sources'] if Path(p).name=='design-si-candidate.obj')
    physical = trimesh.load(physical_path,force='mesh',process=False)
    if not np.array_equal(physical.faces,request.design_faces):
        raise ValueError('PHYSICAL_FACE_IDENTITY_MISMATCH')
    mapping = InProcessCherchiEvaluator(request,spin.SpinEvaluatorConfig(**settings),
        normalization=frozen['normalization'], design_vertices_physical=physical.vertices,
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1', fixed_boundary_reuse=True)
    fixed = json.loads((frozen_path.parent/'fixed.binary64.json').read_text())
    target = (np.asarray(fixed['vertices'])-np.asarray(request.target_origin)) @ np.asarray(request.target_axes).T
    independent_fixed = integrate([[Fraction(float(x)) for x in row] for row in target], fixed['faces'],1.)
    if np.max(np.abs(independent_fixed-mapping.fixed)) > 1e-10:
        raise ValueError('FIXED_MATERIAL_MISMATCH')
    output.mkdir(parents=True,exist_ok=False)
    manifest = dict(samples=old['samples'], layouts=old['layouts'], point_count=81,
        request_sha256=old['request_sha256'], registered_sample_manifest_sha256=hashlib.sha256(old_manifest_path.read_bytes()).hexdigest(),
        settings=settings, normalization=frozen['normalization'], backend='CHERCHI',
        boundary_profile='EXACT_RATIONAL_FIXED_SUBSET_V1',
        scope='ACTUAL_CHERCHI_PRODUCTION_DERIVATIVE_NOT_OPTIMIZATION', case_admission=False)
    (output/'manifest.json').write_text(json.dumps(manifest))
    moments=[]; errors=[]; jacobian=None
    with (output/'events.jsonl').open('x') as events:
        for index,q in enumerate(samples):
            start=time.monotonic()
            try:
                mapped=mapping(q); material=mapping.material(q)
                path=output/f'point-{index:03d}-material.json'
                path.write_text(json.dumps(material,default=_json_value,allow_nan=False))
                source=mapping.last['source']
                branch=dict(surface_type=source['surface_type'], radius=source['truncation_radius'],
                    cap_count=source['cap_face_count'], faces=np.asarray(source['faces']).tolist(),
                    normal_vertices=np.asarray(source['source_normal_vertices']).tolist(),
                    normal_edges=np.asarray(source['source_normal_edges']).tolist(),
                    projective_ray_ids=np.asarray(source['source_projective_rays']).tolist())
                (output/f'point-{index:03d}-source-branch.json').write_text(json.dumps(branch,default=_json_value))
                saved=json.loads(path.read_text())
                snapshot, boundary_mode = rational_solver_candidate_snapshot(
                    saved['variable_rational_vertices'],saved['variable_faces'],
                    fixed_boundary_proof=mapping.fixed_boundary_proof)
                checked_vertices, checked_faces = snapshot.arrays()
                if not np.array_equal(checked_vertices,saved['variable_vertices']) or not np.array_equal(checked_faces,saved['variable_faces']):
                    raise ValueError('BOUNDARY_SNAPSHOT_IDENTITY_MISMATCH')
                exact=[[Fraction(x) for x in row] for row in saved['variable_rational_vertices']]
                independent=integrate(exact,saved['variable_faces'],request.material_density)+independent_fixed
                error=float(np.max(np.abs(independent-mapped['raw_moments'])))
                if error > 1e-10: raise ValueError('INDEPENDENT_MATERIAL_INTEGRATION_FAILED')
                moments.append(independent); errors.append(error)
                if index==0: jacobian=np.asarray(mapped['moment_jacobian'])
                row=dict(index=index,status='RETURNED',seconds=time.monotonic()-start,
                    independent_moments=independent.tolist(), integration_error=error,
                    boundary_validation=boundary_mode)
            except Exception as exc:
                events.write(json.dumps(dict(index=index,status='FAILED',error=f'{type(exc).__name__}: {exc}'))+'\n')
                events.flush(); os.fsync(events.fileno()); raise
            events.write(json.dumps(row)+'\n'); events.flush(); os.fsync(events.fileno())
    checks=[]
    for item in old['layouts']:
        numerical=weights(item['positions'])@np.asarray([moments[i] for i in item['indices']])
        error=float(np.max(np.abs(numerical-jacobian[:,item['direction']])))
        checks.append(dict(**item,error=error,passed=error<=1e-6))
    report=dict(point_count=81,checks=checks,numerical_passed=all(c['passed'] for c in checks),
        independent_error=max(errors),branch_certification=False,case_admission=False)
    (output/'report.json').write_text(json.dumps(report))
    if not report['numerical_passed']: raise ValueError('CHERCHI_FD_FAILED_NO_STEP_RETUNING')


if __name__=='__main__': main(Path(sys.argv[1]))
