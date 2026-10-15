import pytest


def test_seed_protocol_does_not_accept_diagnostic_initialization_or_changed_limits():
    from run_noodles_formal_case import registered_solver
    good=dict(maximum_iterations=60,maximum_restarts=8,maximum_function_evaluations_per_restart=60,
              equation_tolerance=1e-12,deterministic_seed=0,equation_profile='DUAL_MONOMIAL_RAW_V1',
              capture_evaluation_payload=True)
    assert registered_solver(good,1)['deterministic_seed']==1
    assert good['deterministic_seed']==0
    with pytest.raises(ValueError):registered_solver(good,2)
    for key,value in [('maximum_restarts',9),('equation_tolerance',1e-9),
                      ('initial_coefficients',[1.]*10),('capture_evaluation_payload',False)]:
        with pytest.raises(ValueError):registered_solver(dict(good,**{key:value}),0)
