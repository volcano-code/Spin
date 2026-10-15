import importlib
import pytest


def test_offline_classification_requires_own_complete_cherchi_gate():
    module=importlib.import_module('certify_noodles_cherchi_samples')
    report=dict(numerical_passed=True,point_count=81)
    manifest=dict(backend='CHERCHI',point_count=81,request_sha256='expected')
    module.validate_gate(report,manifest,'expected')
    for bad in (dict(report,numerical_passed=False),dict(report,point_count=80)):
        with pytest.raises(ValueError):module.validate_gate(bad,manifest,'expected')
    for bad in (dict(manifest,backend='CGAL'),dict(manifest,request_sha256='other'),dict(manifest,point_count=80)):
        with pytest.raises(ValueError):module.validate_gate(report,bad,'expected')
