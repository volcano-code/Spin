import importlib.util
import unittest


class CorePolicyTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('core_faster_heart'), 'core runner is missing')
        import core_faster_heart
        return core_faster_heart

    def test_core_policy_removes_authority_audits_but_retains_actual_derivative_and_resolution(self):
        original=dict(derivative_method='DISCRETE_FULL_RANK_QUADRIC',icosphere_subdivisions=2,
            maximum_icosphere_subdivisions=2,retain_fixed_normal_provenance=True,
            rational_source_coordinates=True,compact_rational_source_coordinates=True,
            exact_binary64_source_boundary=True,exact_binary64_material_boundary=True,
            cap_radius_policy='GROW_TO_CERTIFIED_V1')
        configured=self.module().core_settings(original)
        self.assertFalse(configured['rational_source_coordinates'])
        self.assertFalse(configured['compact_rational_source_coordinates'])
        self.assertFalse(configured['exact_binary64_material_boundary'])
        self.assertTrue(configured['exact_binary64_source_boundary'])
        self.assertTrue(configured.get('exact_binary64_self_intersection',False),
                        'core geometry must use the predicate confirmed to reject heuristic false positives')
        self.assertEqual(configured['derivative_method'],'DISCRETE_FULL_RANK_QUADRIC')
        self.assertEqual(configured['icosphere_subdivisions'],2)
        self.assertEqual(configured['maximum_icosphere_subdivisions'],2)
        self.assertTrue(original['rational_source_coordinates'])

    def test_another_derivative_method_is_not_silently_relabelled_as_core(self):
        with self.assertRaisesRegex(ValueError,'DISCRETE'):
            self.module().core_settings(dict(derivative_method='LEGACY_CONTINUOUS_IMPLICIT_SHAPE'))


if __name__=='__main__':unittest.main()
