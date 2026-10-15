"""Explicit finite evidence aggregation; no parameter-neighborhood inference."""
from collections import Counter
import math


def summarize(numeric,connections):
    rows=connections['rows'];checks=numeric['checks']
    if (numeric['point_count']!=81 or connections['point_count']!=81 or len(rows)!=81 or
        [r['index'] for r in rows]!=list(range(81))):
        raise ValueError('INCOMPLETE_OR_DUPLICATED_SAMPLE_ROWS')
    expected={(d,h) for d in range(10) for h in (3e-6,1.5e-6,.75e-6)}
    if len(checks)!=30 or {(c['direction'],c['scale']) for c in checks}!=expected:
        raise ValueError('REGISTERED_NUMERIC_LAYOUT_MISMATCH')
    numeric_ok=(numeric['numerical_passed'] is True and
        math.isfinite(numeric['independent_error']) and 0<=numeric['independent_error']<=1e-10 and
        all(c['passed'] is True and math.isfinite(c['error']) and 0<=c['error']<=1e-6 for c in checks))
    equivalent={'IDENTICAL_OR_RENUMBERED','EXACT_COPLANAR_RETRIANGULATION','EXACT_COPLANAR_REFINEMENT'}
    exact_ok=(connections['finite_classification_complete'] is True and
        all(r['status']=='FINITE_SAMPLE_CLASSIFIED' and r['own_construction_replayed'] is True and
            r['same_reference_surface'] is True and r['classification'] in equivalent for r in rows))
    return dict(finite_sample_evidence_passed=numeric_ok and exact_ok,
        numerical_evidence_passed=numeric_ok,finite_oriented_surface_equivalence=exact_ok,
        classification_counts=dict(Counter(r['classification'] for r in rows)),
        maximum_fd_error=max(c['error'] for c in checks),point_count=81,
        neighborhood_certification=False,case_admission=False,
        scope='REGISTERED_FINITE_POINTS_ONLY_NO_RUNTIME_OR_OPTIMIZATION_CLAIM')
