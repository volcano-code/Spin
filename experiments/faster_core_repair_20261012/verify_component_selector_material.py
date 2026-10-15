"""Single registered failure, actual formal Cherchi material after selector repair."""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator,_json_value
from independent_tetra_moments import integrate
HERE=Path(__file__).resolve().parent

def main(output):
    output.mkdir(parents=True,exist_ok=False)
    registration=json.loads((HERE/'noodles-first-core-failures-registration-01/registration.json').read_text())
    config_path=Path(registration['mapping_config']);raw=config_path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=registration['mapping_config_sha256']:raise ValueError('CHANGED_CONFIG')
    config=json.loads(raw);payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**config['evaluator'])
    q=np.asarray(registration['failures']['source_side']['record']['coefficients'])
    adapter=InProcessCherchiEvaluator(request,settings);start=time.monotonic()
    mapped=adapter(q);material=adapter.material(q)
    (output/'material.json').write_text(json.dumps(material,default=_json_value,allow_nan=False))
    original=json.loads((HERE/'core-registered-failure-capture-job-01/output/source_side.json').read_text())['snapshots'][0]
    source=adapter.last['source']
    same_vertices=np.array_equal(source['vertices'],original['source_vertices'])
    same_faces=np.array_equal(source['faces'],original['source_faces'])
    independent=integrate(adapter.last['exact_vertices'],adapter.last['faces'],request.material_density)
    error=float(np.max(np.abs(independent-mapped['variable_raw_moments'])))
    report=dict(same_source_vertices=bool(same_vertices),same_source_faces=bool(same_faces),
        independent_moment_error=error,derivative_defined=mapped['derivative_defined'],
        seconds=time.monotonic()-start,call_counts=adapter.call_counts,core_admission=False,
        scope='Registered state actual material; changed selector derivative validation still required')
    (output/'report.json').write_text(json.dumps(report,allow_nan=False))

if __name__=='__main__':main(Path(sys.argv[1]))
