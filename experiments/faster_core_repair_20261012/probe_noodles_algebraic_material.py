"""One actual CGAL material map after registered source arithmetic correction."""
from dataclasses import replace
import hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'faster_heart_contact_20261011'))
from make_it_stand.spin import SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig
from rounded_solver_journal import run_durable_fullrank,save
report_path=HERE/'noodles-interval-algebraic-source-derivative-job-01/output/report.json'
gate=json.loads(report_path.read_text())
if not gate['source_checks_passed'] or gate['point_count']!=81:raise ValueError('SOURCE_GATE_REQUIRED')
frozen=json.loads((HERE/'noodles-balanced-request-freeze-job-01/output/request.json').read_text())
row=next(x for x in frozen['requests'] if x['settings']['deterministic_seed']==0)
settings=SpinSolverConfig(**row['settings']);payload=dict(row)
payload['settings']=replace(settings,maximum_iterations=1,maximum_restarts=1,maximum_function_evaluations_per_restart=1)
evaluator=SpinEvaluatorConfig(**frozen['evaluator_settings']);evaluator=replace(evaluator,algebraic_finite_map=True)
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
save(out/'purpose.json',dict(role='ONE_SEED_INITIAL_ACTUAL_MATERIAL_DIAGNOSIS',full_solver_limits=row['settings'],diagnostic_solver_limits=payload['settings'].__dict__,explicit_profile='PAPER_INTERVAL_ALGEBRAIC_FINITE_CORE_DIAGNOSIS',source_gate_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),optimization_admission=False))
result=run_durable_fullrank(SpinItFasterRequest(**payload),evaluator,out/'run')
save(out/'report.json',dict(role='ONE_SEED_INITIAL_ACTUAL_MATERIAL_DIAGNOSIS',solver_status=result['status'],function_evaluations=result['function_evaluations'],optimization_admission=False,production_derivative_batch_admission=False))
