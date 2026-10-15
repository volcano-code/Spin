from copy import deepcopy
import pytest


def test_only_certified_wall_metadata_and_roundoff_fixed_moments_may_differ():
    from noodles_formal_case import compare_requests
    old=dict(minimum_wall_thickness=0.,fixed_shell_moments=[1.]*10,
             center_of_mass_weight=1.,design_vertices=[[1.,2.,3.]],settings={'seed':0})
    new=deepcopy(old);new['minimum_wall_thickness']=.02;new['fixed_shell_moments'][0]+=1e-15
    assert compare_requests(old,new,.02)['equivalent']
    for key,value in [('center_of_mass_weight',2.),('design_vertices',[[1.,2.,3.1]]),
                      ('minimum_wall_thickness',.03),('fixed_shell_moments',[1.01]*10)]:
        bad=deepcopy(new);bad[key]=value
        with pytest.raises(ValueError):compare_requests(old,bad,.02)
