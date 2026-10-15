"""Locked small top: production Cherchi seed-only solve, seal and fresh reload."""
from dataclasses import asdict
import json
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from make_it_stand.verification.spin_bundle import save_spin_artifact
from make_it_stand.verification.cherchi_inprocess import solve_cherchi_inprocess_and_seal_faster

HERE=Path(__file__).resolve().parent
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
config=json.loads((HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())
frozen=json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())
r=config['request'];fixed=frozen['fixed_material']
problem=SpinArtifactProblem(original=(r['outer_vertices'],r['outer_faces']),
    fixed_regions=[('shell',(fixed['vertices'],fixed['faces']),fixed['density'])],
    design_regions=[('design',(r['design_vertices'],r['design_faces']),r['material_density'])],
    length_unit='m',mass_unit='kg',minimum_outer_wall=0.05,
    target=dict(mode='SPINNING_TOP',axis=[0.,0.,1.],origin=r['target_origin'],support_points=[]))
settings=spin.SpinEvaluatorConfig(**config['evaluator'])
solver=spin.SpinSolverConfig(**r['settings'])
with (out/'registered-configuration.json').open('x') as stream:
    json.dump(dict(problem_digest=problem.input_digest,normalization=problem.normalization,
        evaluator=asdict(settings),solver=asdict(solver),weight=1.,
        role='FORMAL_CHERCHI_NORMALIZED_SMALL_TOP_CLOSURE_ONLY',
        initialization='REQUEST_SEED_ONLY',complete_core_admission=False),stream,sort_keys=True)
start=time.monotonic()
result=solve_cherchi_inprocess_and_seal_faster(problem,evaluator_settings=settings,
    solver_settings=solver,output_directory=out/'run',center_of_mass_weight=1.)
if result.status!='CONVERGED' or result.artifact is None:
    with (out/'report.json').open('x') as stream:
        json.dump(dict(status=result.status,seal_error=result.seal_error,
            closure_passed=False,seconds=time.monotonic()-start),stream)
    raise SystemExit(1)
run=result.identity['native_run']
if not any(row.get('accepted_geodesic_step',0)>0 for row in run['trace']):
    raise ValueError('NO_ACCEPTED_UPDATE')
if run['raw_equation_residual']>1e-12:raise ValueError('RAW_RESIDUAL_GATE')
saved=save_spin_artifact(result.artifact,out/'material')
# save_spin_artifact already returns the independently reloaded reader object.
loaded=saved
with (out/'report.json').open('x') as stream:
    json.dump(dict(status='FORMAL_CHERCHI_SMALL_CLOSURE_PASSED',closure_passed=True,
        raw_equation_residual=run['raw_equation_residual'],seconds=time.monotonic()-start,
        scientific_digest=result.scientific_digest,material_digest=loaded.scientific_digest,
        source='NATIVE_SEED_INITIALIZATION',core_admission=False,
        scope='Normalized small formal closure; derivative/repeat/performance gates pending'),stream,sort_keys=True)
print('FORMAL_CHERCHI_SMALL_CLOSURE_PASSED',flush=True)
