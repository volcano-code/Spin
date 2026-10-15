"""Explicitly scoped equivalence check; preserves every geometric fact."""
def compare_boundary_facts(baseline, candidate):
    if candidate.get('backend') != baseline.get('backend','') + '/conservative-exact-diagonal-slabs-v1':
        raise ValueError('UNREGISTERED_BOUNDARY_BACKEND_CHANGE')
    old, new = dict(baseline), dict(candidate)
    old.pop('backend'); new.pop('backend')
    old_count, new_count = old.pop('exact_pair_tests'), new.pop('exact_pair_tests')
    if old != new:
        raise ValueError('BOUNDARY_GEOMETRIC_FACTS_CHANGED')
    if not 0 <= new_count <= old_count:
        raise ValueError('PREDICATE_COUNT_NOT_REDUCED')
    return dict(geometric_facts_identical=True, baseline_exact_pair_tests=old_count,
                candidate_exact_pair_tests=new_count, saved_exact_pair_tests=old_count-new_count)
