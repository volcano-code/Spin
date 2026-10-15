from copy import deepcopy
import pytest


def test_only_complete_numeric_and_exact_equivalent_samples_pass():
    from finite_sample_summary import summarize
    numeric=dict(point_count=81,numerical_passed=True,independent_error=1e-17,
        checks=[dict(direction=d,scale=h,error=1e-10,passed=True)
                for d in range(10) for h in (3e-6,1.5e-6,.75e-6)])
    connections=dict(point_count=81,finite_classification_complete=True,rows=[
        dict(index=i,status='FINITE_SAMPLE_CLASSIFIED',classification='EXACT_COPLANAR_RETRIANGULATION',
             own_construction_replayed=True,same_reference_surface=True) for i in range(81)])
    result=summarize(numeric,connections)
    assert result['finite_sample_evidence_passed']
    assert not result['neighborhood_certification'] and not result['case_admission']
    for field,value in [('classification','REAL_CONSTRUCTION_BRANCH_CHANGE'),
                        ('own_construction_replayed',False),('same_reference_surface',False)]:
        bad=deepcopy(connections);bad['rows'][1][field]=value
        assert not summarize(numeric,bad)['finite_sample_evidence_passed']
    bad=deepcopy(numeric);bad['checks'][0]['error']=2e-6
    assert not summarize(bad,connections)['finite_sample_evidence_passed']
    bad=deepcopy(connections);bad['rows'][1]['index']=0
    with pytest.raises(ValueError):summarize(numeric,bad)
