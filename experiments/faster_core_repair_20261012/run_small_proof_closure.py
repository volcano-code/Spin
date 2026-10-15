"""Seeded CGAL integration regression with same-call terminal/bundle reload.

Not a Heart/Noodles result or general backend admission. Reuses existing mesh
validation and bundle I/O; does not re-evaluate the accepted terminal quadric.
"""
from dataclasses import asdict
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.spin_problem import MeshSnapshot, partition
from make_it_stand.verification.spin_bundle import _Writer, _Reader, SCHEMA
from independent_tetra_moments import integrate
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'faster_heart_contact_20261011'))
from rounded_solver_journal import run_durable_fullrank, save

out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
config=json.loads((HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())
frozen=json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())
payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
request=spin.SpinItFasterRequest(**payload)
settings=spin.SpinEvaluatorConfig(**dict(config['evaluator'],exact_binary64_material_boundary=True))
save(out/'configuration.json',dict(request=asdict(request),evaluator=asdict(settings),
    role='SMALL_FORMAL_INTERFACE_SAME_CALL_TERMINAL_RELOAD_REGRESSION',initialization='FROZEN_REQUEST_SEED_ONLY'))
start=time.monotonic()
result=run_durable_fullrank(request,settings,out/'run',capture_material_proof=True)
solve_seconds=time.monotonic()-start
if result['status']!='CONVERGED':
    save(out/'report.json',dict(status=result['status'],closure_passed=False,solve_seconds=solve_seconds,optimization_case_admission=False))
    sys.exit(1)
ref=json.loads((out/'run/terminal-material-reference.json').read_text())
terminal=json.loads((out/'run'/ref['evaluation_file']).read_text())['evaluation']
proof=json.loads(terminal['material_proof_json'])
def rational_rows(points):
    return [[f'{x.numerator}/{x.denominator}' for x in (Fraction(v) for v in row)] for row in points]
def to_target(points):return (np.asarray(points)-request.target_origin)@np.asarray(request.target_axes)
retained=MeshSnapshot.from_rationals(rational_rows(proof['exact_material_vertices']),terminal['material_faces'])
design=MeshSnapshot((to_target(request.design_vertices),request.design_faces))
outer=MeshSnapshot((to_target(request.outer_vertices),request.outer_faces))
fixed_config=frozen['fixed_material']
fixed=MeshSnapshot((to_target(fixed_config['vertices']),fixed_config['faces']))
native=spin._native_module()
cavity_packet=native.libigl_exact_boolean_construction_from_rationals(
    rational_rows(to_target(request.design_vertices)),np.asarray(request.design_faces),
    rational_rows(proof['exact_material_vertices']),np.asarray(terminal['material_faces']),'difference')
if not cavity_packet['succeeded']:raise ValueError('COMPLEMENT_CONSTRUCTION_FAILED')
cavity=MeshSnapshot.from_rationals(rational_rows(cavity_packet['vertices']),cavity_packet['faces'])
materials=dict(outer=outer,design=design,fixed=fixed,retained=retained,cavity=cavity)
bundle=out/'material';bundle.mkdir(exist_ok=False)
fd=os.open(bundle,os.O_RDONLY|os.O_DIRECTORY)
try:
    writer=_Writer(fd)
    records={name:writer.mesh(mesh) for name,mesh in materials.items()}
    writer.write('manifest.json',json.dumps(dict(schema=SCHEMA,roles=records,terminal_reference=ref),sort_keys=True).encode())
    os.fsync(fd)
finally:os.close(fd)
# A fresh Reader reconstructs the saved exact regions; no writer object cache.
fd=os.open(bundle,os.O_RDONLY|os.O_DIRECTORY)
try:
    reader=_Reader(fd)
    manifest=json.loads(reader.read('manifest.json'))
    loaded={name:reader.mesh(record,schema=manifest['schema']) for name,record in manifest['roles'].items()}
finally:os.close(fd)
design_partition=partition(loaded['design'].exact_region(),[loaded['retained'].exact_region(),loaded['cavity'].exact_region()])
outer_partition=partition(loaded['outer'].exact_region(),[loaded['fixed'].exact_region(),loaded['retained'].exact_region(),loaded['cavity'].exact_region()])
def independent(mesh,density):
    vertices,faces=mesh.arrays()
    exact=mesh.rational_vertices
    if exact is not None:vertices=[[Fraction(x) for x in p] for p in exact]
    return integrate(vertices,faces,density)
variable=independent(loaded['retained'],request.material_density)
total=variable+independent(loaded['fixed'],fixed_config['density'])
error=float(np.max(np.abs(total-np.asarray(terminal['total_raw_moments']))))
if error>1e-6 or result['raw_equation_residual']>1e-12:raise ValueError('TERMINAL_NUMERICAL_CHECK_FAILED')
if not any(row.get('accepted_geodesic_step',0)>0 for row in result['trace']):raise ValueError('NO_ACCEPTED_OPTIMIZATION_UPDATE')
save(out/'report.json',dict(status='SMALL_SAME_CALL_CLOSURE_PASSED',closure_passed=True,
    solve_seconds=solve_seconds,total_seconds=time.monotonic()-start,
    independent_total_moments=total.tolist(),maximum_moment_error=error,
    raw_equation_residual=result['raw_equation_residual'],design_partition=design_partition,
    outer_partition=outer_partition,terminal_quadric_reevaluations=0,
    optimization_case_admission=False,scope='Small seeded CGAL integration regression only; no case repeat or physics admission'))
print('SMALL_SAME_CALL_CLOSURE_PASSED',flush=True)
