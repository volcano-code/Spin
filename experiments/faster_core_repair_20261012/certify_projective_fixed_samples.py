import json
from pathlib import Path
import sys
from certify_formal_fixed_samples import main

HERE = Path(__file__).resolve().parent
inputs = HERE / 'projective-formal-pipeline-job-01/output'
registered = json.loads((inputs / 'registered-inputs.json').read_text())
if registered['settings']['paper_projective_hyperboloid'] is not True:
    raise ValueError('PROJECTIVE_STRATEGY_NOT_REGISTERED')
main(Path(sys.argv[1]), inputs / 'derivatives', registered['request']['design_vertices'])
