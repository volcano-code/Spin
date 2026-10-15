"""Original registered state, unchanged q and three original scales."""
import hashlib,json,sys
from pathlib import Path
from probe_formal_fixed_derivatives import main
HERE=Path(__file__).resolve().parent
registration=json.loads((HERE/'noodles-first-core-failures-registration-01/registration.json').read_text())
raw=Path(registration['mapping_config']).read_bytes()
if hashlib.sha256(raw).hexdigest()!=registration['mapping_config_sha256']:raise ValueError('CHANGED_CONFIG')
config=json.loads(raw)
main(Path(sys.argv[1]),dict(request=config['request'],settings=config['evaluator'],
    states=dict(source_side=registration['failures']['source_side']['record']['coefficients'])),normalize=False)
