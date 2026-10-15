"""Independent hand-derived fixtures; no native production derivative imports."""
import importlib.util
from pathlib import Path
import sys
import mpmath as mp
import pytest

HERE=Path(__file__).resolve().parent
def reference():
    path=HERE/'independent_directional_reference.py'
    assert path.exists(), 'independent high-precision reference missing'
    spec=importlib.util.spec_from_file_location('independent_directional_reference',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module

def test_sphere_normal_motion_matches_radius_derivative():
    r=reference()
    with mp.workdps(70):
        q=list(map(mp.mpf,[-1,0,0,0,1,1,1,0,0,0]))
        x=r.finite_point(q,[1,0,0])
        assert max(abs(a-b) for a,b in zip(x,[1,0,0]))<mp.mpf('1e-60')
        dx=mp.diff(lambda t:r.finite_point([q[0]+t,*q[1:]],[1,0,0])[0],0)
        assert abs(dx+mp.mpf('0.5'))<mp.mpf('1e-60')

def test_ellipsoid_axis_motion_is_not_sphere_motion():
    r=reference()
    with mp.workdps(70):
        q=[mp.mpf(-1),0,0,0,mp.mpf(1)/4,mp.mpf(1)/9,mp.mpf(1)/16,0,0,0]
        assert abs(r.finite_point(q,[1,0,0])[0]-2)<mp.mpf('1e-60')
        assert abs(mp.diff(lambda t:r.finite_point([q[0]+t,*q[1:]],[1,0,0])[0],0)+1)<mp.mpf('1e-60')

def test_line_plane_intersection_motion():
    r=reference();D=r.Dual
    point=r.line_plane([[D(0),D(0),D(-1)],[D(0),D(0),D(1)],
                        [D(0),D(0),D(0,1)],[D(1),D(0),D(0,1)],[D(0),D(1),D(0,1)]])
    assert point[2].v==0 and point[2].d==1

def test_signed_tetrahedron_mass_motion():
    r=reference();D=r.Dual
    # Axis tetrahedron volume a/6, x first moment a^2/24, x^2 moment a^3/60.
    vertices=[[D(0),D(0),D(0)],[D(2,1),D(0),D(0)],[D(0),D(1),D(0)],[D(0),D(0),D(1)]]
    faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]]
    with mp.workdps(70):
        result=r.tetra_moments(vertices,faces,1)
        for i,derivative in [(0,mp.mpf(1)/6),(1,mp.mpf(1)/6),(4,mp.mpf(1)/5)]:
            assert abs(result[i].d-derivative)<mp.mpf('1e-60')

def test_singular_line_plane_is_rejected():
    r=reference();D=r.Dual
    with pytest.raises(ValueError,match='SINGULAR'):
        r.line_plane([[D(0),D(0),D(1)],[D(1),D(0),D(1)],
                      [D(0),D(0),D(0)],[D(1),D(0),D(0)],[D(0),D(1),D(0)]])

def test_hyperboloid_truncation_root_implicit_derivative():
    r=reference()
    with mp.workdps(70):
        q=list(map(mp.mpf,[-1,0,0,0,1,1,-1,0,0,0]))
        p=[mp.sqrt(3),mp.mpf(0),mp.sqrt(2)]
        point,J,condition=r.root_point_and_jacobian(q,[0,1,0],mp.sqrt(5),p)
        assert abs(J[0,0]+1/(4*mp.sqrt(3)))<mp.mpf('1e-60')
        assert abs(J[2,0]-1/(4*mp.sqrt(2)))<mp.mpf('1e-60')
        assert condition<100

def test_wrong_jacobian_and_uncertain_reference_cannot_pass():
    r=reference()
    assert r.classify([[1.001]],[[1]],[[1]],True)=='INCONSISTENT'
    assert r.classify([[1]],[[1]],[[1]],True)=='SUPPORTED'
    assert r.classify([[1]],[[1]],[[1.01]],True)=='INCONCLUSIVE'
    assert r.classify([[1]],[[1]],[[1]],False)=='INCONCLUSIVE'
