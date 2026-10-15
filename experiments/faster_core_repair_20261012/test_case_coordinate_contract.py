"""Length conversion must preserve the frozen Faster objective, not retune it."""
import importlib
import numpy as np
import pytest


def objective(r, weight):
    m,x,y,z,xx,yy,zz,yz,xz,xy = r
    a=yy+zz-(y*y+z*z)/m
    b=xx+zz-(x*x+z*z)/m
    c=xx+yy-(x*x+y*y)/m
    d=-xy+x*y/m
    return (a*a+b*b+2*d*d)/(c*c)+weight*z/m


@pytest.mark.parametrize('length,density', [(2.57,1.), (.025,1200.), (1.,1.)])
def test_frozen_target_frame_objective_preserved(length,density):
    contract=importlib.import_module('case_coordinate_contract')
    model=np.array([3., .2, -.3, .7, 1.2, 1.5, 2., .03, -.06, .04])
    original=model.copy()
    normalized=contract.normalize_target_moments(model,length,density)
    weight=contract.normalized_faster_weight(1.,length)
    assert objective(normalized,weight)==pytest.approx(objective(model,1.),abs=2e-14)
    np.testing.assert_array_equal(model,original)
    # A bare label change or unchanged numeric weight changes the problem.
    if length != 1.:
        assert abs(objective(normalized,1.)-objective(model,1.))>1e-3


@pytest.mark.parametrize('length', [0., -1., float('nan'), float('inf')])
def test_unrepresentable_coordinate_contract_is_rejected(length):
    contract=importlib.import_module('case_coordinate_contract')
    with pytest.raises(ValueError):
        contract.normalized_faster_weight(1.,length)


def test_wrong_moment_shape_and_nonfinite_are_rejected():
    contract=importlib.import_module('case_coordinate_contract')
    for values in ([1.]*9,[1.]*9+[float('nan')]):
        with pytest.raises(ValueError):
            contract.normalize_target_moments(values,2.57,1.)


def test_quadric_unit_change_preserves_material_side_and_chain_rule():
    contract=importlib.import_module('case_coordinate_contract')
    q=np.array([-.4,.3,.1,.5,.07,.44,-.39,-.07,.18,-.24])
    q/=np.linalg.norm(q)
    length=2.57
    p,jac=contract.normalize_target_quadric(q,length)
    def polynomial(c,x):
        a,b,z=x
        return np.dot(c,[1,a,b,z,a*a,b*b,z*z,2*b*z,2*a*z,2*a*b])
    divisor=np.linalg.norm(q*np.array([1,length,length,length,*([length**2]*6)]))
    for x in ([.2,-.3,.7],[-1.,.8,.1],[0,0,0]):
        assert polynomial(p,np.array(x)/length)==pytest.approx(polynomial(q,x)/divisor,abs=1e-14)
    for direction in np.eye(10):
        h=1e-6
        plus=contract.normalize_target_quadric(q+h*direction,length)[0]
        minus=contract.normalize_target_quadric(q-h*direction,length)[0]
        np.testing.assert_allclose((plus-minus)/(2*h),jac@direction,rtol=1e-7,atol=1e-9)
    np.testing.assert_allclose(jac@q,0,atol=1e-15)
    np.testing.assert_allclose(p@jac,0,atol=1e-15)
