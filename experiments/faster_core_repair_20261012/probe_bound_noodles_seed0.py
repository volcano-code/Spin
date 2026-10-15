"""Exactly one seed-initial production evaluation; NOT a complete solve."""
from dataclasses import replace
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'faster_heart_contact_20261011'))
from make_it_stand.spin import SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig
from rounded_solver_journal import run_durable_fullrank,save
source=Path(sys.argv[1]);output=Path(sys.argv[2]);output.mkdir(parents=True,exist_ok=False)
frozen=json.loads(source.read_text());row=next(r for r in frozen['requests'] if r['settings']['deterministic_seed']==0)
settings=SpinSolverConfig(**row['settings'])
# Sampling uses only mode/seed; limits bound this diagnosis, never real optimization.
payload=dict(row);payload['settings']=replace(settings,maximum_iterations=1,
    maximum_restarts=1,maximum_function_evaluations_per_restart=1)
request=SpinItFasterRequest(**payload)
save(output/'purpose.json',dict(role='ONE_SEED_INITIAL_PRODUCTION_MAP_DIAGNOSIS_ONLY',
    production_request=str(source),full_solver_limits=row['settings'],
    diagnostic_solver_limits=payload['settings'].__dict__,optimization_admission=False))
result=run_durable_fullrank(request,SpinEvaluatorConfig(**frozen['evaluator_settings']),output/'run')
save(output/'report.json',dict(role='ONE_SEED_INITIAL_PRODUCTION_MAP_DIAGNOSIS_ONLY',
    solver_status=result['status'],function_evaluations=result['function_evaluations'],
    optimization_admission=False,derivative_batch_admission=False))
