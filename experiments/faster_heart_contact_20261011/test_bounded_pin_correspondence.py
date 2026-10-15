import importlib.util
import unittest


class WholeBoundaryTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('bounded_pin_correspondence'))
        import bounded_pin_correspondence
        return bounded_pin_correspondence

    def triangle(self,z='0/1'):
        return dict(vertices=[['0/1','0/1',z],['1/1','0/1',z],['0/1','1/1',z]],faces=[[0,1,2]])

    def test_full_parallel_triangle_distance_is_certified_not_sampled(self):
        result=self.module().boundary_distance_certificate(self.triangle(),self.triangle('1/100'),'1/20')
        self.assertEqual(result['status'],'CERTIFIED')
        self.assertEqual(result['bidirectional_distance_squared_upper_bound_exact'],'1/10000')
        self.assertEqual(result['coverage_scope'],'ALL_TRIANGLE_POINTS_VIA_EXACT_CONVEX_COVER')

    def test_remeshed_face_cover_is_allowed_only_when_whole_triangle_covered(self):
        target=dict(vertices=[['0/1']*3,['1/1','0/1','0/1'],['0/1','1/1','0/1'],['1/2','1/2','0/1']],
                    faces=[[0,1,3],[0,3,2]])
        result=self.module().boundary_distance_certificate(self.triangle(),target,'1/1000')
        self.assertEqual(result['status'],'CERTIFIED')
        target['faces']=[[0,1,3]]
        result=self.module().boundary_distance_certificate(self.triangle(),target,'1/1000')
        self.assertEqual(result['status'],'UNKNOWN')
        self.assertFalse(result['geometry_admission_granted'])

    def test_too_distant_or_projection_degenerate_is_unknown_not_pass(self):
        result=self.module().boundary_distance_certificate(self.triangle(),self.triangle('1/1'),'1/100')
        self.assertEqual(result['status'],'UNKNOWN')
        perpendicular=dict(vertices=[['0/1']*3,['0/1','1/1','0/1'],['0/1','0/1','1/1']],faces=[[0,1,2]])
        result=self.module().boundary_distance_certificate(self.triangle(),perpendicular,'1/100')
        self.assertEqual(result['status'],'UNKNOWN')

    def test_zero_moment_error_cannot_replace_full_boundary_certificate(self):
        module=self.module()
        self.assertTrue(callable(getattr(module,'numerical_pin_gate',None)))
        report=module.numerical_pin_gate(dict(status='UNKNOWN'),['1/1']+['0/1']*9,
                                        ['1/1']+['0/1']*9,'0/1','1/1')
        self.assertEqual(report['status'],'UNKNOWN')
        self.assertFalse(report['geometry_admission_granted'])

    def test_wrong_difference_volume_or_moment_error_is_rejected(self):
        module=self.module()
        self.assertTrue(callable(getattr(module,'numerical_pin_gate',None)))
        boundary=dict(status='CERTIFIED',bidirectional_distance_squared_upper_bound_exact='0/1')
        with self.assertRaisesRegex(ValueError,'DIFFERENCE_VOLUME'):
            module.numerical_pin_gate(boundary,['1/1']+['0/1']*9,['1/1']+['0/1']*9,'1/10','1/1')
        report=module.numerical_pin_gate(boundary,['1/1']+['0/1']*9,['999/1000']+['0/1']*9,'1/1000','1/1')
        self.assertEqual(report['status'],'REJECTED')

if __name__=='__main__':unittest.main()
