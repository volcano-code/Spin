"""Saved registered projective materials only; no new Boolean or source calls."""
import hashlib,json,sys
from pathlib import Path
from certify_formal_fixed_samples import main
HERE=Path(__file__).resolve().parent
registration=json.loads((HERE/'noodles-first-core-failures-registration-01/registration.json').read_text())
raw=Path(registration['mapping_config']).read_bytes()
if hashlib.sha256(raw).hexdigest()!=registration['mapping_config_sha256']:raise ValueError('CHANGED_CONFIG')
config=json.loads(raw)
main(Path(sys.argv[1]),HERE/'component-selector-derivatives-job-01/output',
    config['request']['design_vertices'],('source_side',))
