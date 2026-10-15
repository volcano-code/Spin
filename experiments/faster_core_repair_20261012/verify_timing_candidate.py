"""Verify the actual loaded candidate's outputs and stage propagation."""
import json
from pathlib import Path
import numpy as np


def verify(candidate,baseline):
    report=json.loads((candidate/'report.json').read_text())
    if any('exception' in row for row in report['rows']):raise ValueError('CANDIDATE_DIAGNOSTIC_FAILURE')
    checks=[]
    for state in ('ellipsoid','one_sheet','two_sheet'):
        before=json.loads((baseline/(state+'-fullrank-evaluation.json')).read_text())
        after=json.loads((candidate/(state+'-fullrank-evaluation.json')).read_text())
        if before != after:raise ValueError('INSTRUMENTATION_CHANGED_SCIENTIFIC_OUTPUT:'+state)
        timing=json.loads((candidate/(state+'-cgal-0-record.json')).read_text())['fullrank_phases']
        if timing['incomplete_stages']:raise ValueError('INCOMPLETE_NATIVE_TIMING:'+state)
        for name in ('canonical-source-motion','canonical-material-transport'):
            if timing['stages'].get(name,{}).get('calls')!=1:raise ValueError('MISSING_CANONICAL_STAGE:'+name)
        comparison=json.loads((candidate/(state+'-cherchi-0-record.json')).read_text())
        if comparison['raw_moment_max_error']>1e-10 or comparison['jacobian_max_error']>1e-6:
            raise ValueError('BACKEND_COMPONENT_MISMATCH:'+state)
        checks.append(dict(state=state,scientific_output_identical=True,canonical_timing_complete=True))
    return checks


if __name__=='__main__':
    import sys
    print(json.dumps(verify(Path(sys.argv[1]),Path(sys.argv[2]))))
