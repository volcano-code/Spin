import pytest
from boundary_candidate_comparison import compare_boundary_facts


def test_only_registered_backend_and_pair_count_change_is_allowed():
    old = dict(backend='exact', exact_pair_tests=20, boundary_valid=True, input_digest='abc')
    new = dict(old, backend='exact/conservative-exact-diagonal-slabs-v1', exact_pair_tests=12)
    result = compare_boundary_facts(old,new)
    assert result['saved_exact_pair_tests'] == 8
    for key,value in [('boundary_valid',False),('input_digest','other'),('extra',1),
                      ('backend','unregistered')]:
        with pytest.raises(ValueError):
            compare_boundary_facts(old,dict(new,**{key:value}))


def test_increased_predicate_work_is_not_a_speedup_certificate():
    old = dict(backend='exact', exact_pair_tests=20, boundary_valid=True)
    with pytest.raises(ValueError):
        compare_boundary_facts(old,dict(old,backend='exact/conservative-exact-diagonal-slabs-v1',exact_pair_tests=21))
