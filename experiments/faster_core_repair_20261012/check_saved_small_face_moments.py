"""Independent exact integration of the saved same-call repaired material."""
import hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
import numpy as np
from independent_tetra_moments import integrate
HERE=Path(__file__).resolve().parent
path=HERE/'core-small-face-verification-job-01/output/actual-evaluation.json'
saved=json.loads(path.read_text())
config_path=HERE/'noodles-registered-core-seed-0-job-01/output/run/mapping-config.json'
rho=json.loads(config_path.read_text())['request']['material_density']
start=time.monotonic()
binary=integrate(saved['material_vertices'],saved['material_faces'],rho)
exact=integrate([[Fraction(x) for x in row] for row in saved['material_proof_attachment']['exact_material_vertices']],
                saved['material_faces'],rho)
expected=np.asarray(saved['variable_raw_moments'])
report=dict(schema='independent-saved-fixed-state-moments-v1',production_admission=False,optimization_invocations=0,
            material_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
            binary64_exact_integration_error=float(max(abs(binary-expected))),
            rational_exact_integration_error=float(max(abs(exact-expected))),
            derivative_defined=saved['derivative_defined'],
            production_derivative_certified=False,seconds=time.monotonic()-start)
report['independent_variable_moments_passed']=max(report['binary64_exact_integration_error'],
                                              report['rational_exact_integration_error'])<=1e-10
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
with (out/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report),flush=True)
if not report['independent_variable_moments_passed']:raise ValueError('INDEPENDENT_VARIABLE_MOMENTS_FAILED')
