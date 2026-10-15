import importlib
import mpmath as mp
import pytest


def test_ideal_edge_keeps_both_roots_in_parameter_order():
    module=importlib.import_module('independent_projective_reference')
    with mp.workdps(70):
        q=[-1,0,0,0,1,1,-1,0,0,0]
        roots=module.ideal_edge(q,[-2,0,1],[2,0,1])
        assert roots==[(mp.mpf('.25'),[-1,0,-1]),(mp.mpf('.75'),[1,0,-1])]
        with pytest.raises(ValueError):
            module.ideal_edge(q,[1,0,1],[2,0,1])


def test_ray_hits_sphere_not_analytic_quadric_and_has_correct_radius_motion():
    module=importlib.import_module('independent_projective_reference')
    with mp.workdps(70):
        point=module.truncate_ray([1,0,0],[mp.mpf('.5'),0,mp.mpf('-.5')],3)
        root=mp.sqrt(17)
        expected=[(1+root)/2,0,(1-root)/2]
        assert max(abs(a-b) for a,b in zip(point,expected))<mp.mpf('1e-65')
        derivative=[mp.diff(lambda r:module.truncate_ray([1,0,0],[1,0,-1],r)[i],3) for i in range(3)]
        assert max(abs(a-b) for a,b in zip(derivative,[3/root,0,-3/root]))<mp.mpf('1e-65')
        with pytest.raises(ValueError):module.truncate_ray([3,0,0],[1,0,0],3)


def test_remote_cap_preserves_tangent_coordinates_and_rejects_ambiguous_side():
    module=importlib.import_module('independent_projective_reference')
    with mp.workdps(70):
        loop=[[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]]
        expected=[[-1,-1,6],[1,-1,6],[1,1,6],[-1,1,6]]
        assert module.remote_cap(loop,3)==expected
        assert module.remote_cap(list(reversed(loop)),3)==list(reversed(expected))
        motion=[mp.diff(lambda r:module.remote_cap(loop,r)[0][i],3) for i in range(3)]
        assert motion==[0,0,2]
        with pytest.raises(ValueError):module.remote_cap([[x,y,0] for x,y,z in loop],3)
        with pytest.raises(ValueError):module.remote_cap([[0,0,1],[1,0,1],[2,0,1]],3)


def test_projective_recipe_replays_moving_normal_and_both_ray_branches():
    module=importlib.import_module('independent_projective_reference')
    with mp.workdps(70):
        points=module.projective_points([-1,0,0,0,-1,1,1,0,0,0],
            [[-2,0,1],[2,0,1]], [2,-1,-1],
            [[-3,0,1,0],[2,0,1,0],[2,0,1,1]],3,[],3)
        assert points[0]==[0,0,1]
        assert points[1][0]>0 and points[2][0]<0
        assert points[1][2]==points[2][2]
        for point in points[1:]:
            assert abs(sum(x*x for x in point)-9)<mp.mpf('1e-65')
        with pytest.raises(ValueError):
            module.projective_points([-1,0,0,0,-1,1,1,0,0,0],
                [[-2,0,1],[2,0,1]],[-1],[[2,0,1,0]],1,[],3)


def test_projective_recipe_assembles_registered_extruded_cap():
    module=importlib.import_module('independent_projective_reference')
    with mp.workdps(70):
        normals=[[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]]
        rays=[[-1]*4 for _ in normals]+[[-2,i,0,i] for i in range(4)]
        points=module.projective_points([-1,0,0,0,1,1,1,0,0,0],
            normals,[0,1,2,3,-1,-1,-1,-1],rays,4,[],3)
        assert len(points)==8
        for before,after in zip(points[:4],points[4:]):
            assert before[:2]==after[:2] and after[2]==6
