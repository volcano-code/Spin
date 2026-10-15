"""The external geometry entry must select the original native objective."""
import numpy as np
import pytest
from make_it_stand import spin


def test_original_objective_mode_mismatch_rejected_before_geometry():
    native = spin._native_module()
    settings = getattr(native, 'Spin2014QuadricObjectiveSettings', None)
    assert settings is not None, 'formal unified build must enable the original objectives'
    objective = settings()
    objective.objective = native.Spin2014Objective.YOYO
    request = native.SpinDesignRequest()
    request.mode = native.SpinMode.SPINNING_TOP
    request.settings.equation_profile = native.SpinEquationProfile.DUAL_MONOMIAL_RAW_V1
    def geometry(q):
        raise AssertionError('mode mismatch must fail before geometry')
    with pytest.raises(ValueError, match='mode'):
        native.solve_spin_it_faster_external_moment_evaluator(request, geometry,
            np.empty((0, 10)), objective_settings=objective)
