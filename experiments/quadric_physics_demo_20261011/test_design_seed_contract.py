import unittest
import core_boolean_mesh as core


class DesignSeedContractTests(unittest.TestCase):
    def select(self,expression,parts):
        fn=getattr(core,'manifest_design_seed',None)
        self.assertTrue(callable(fn),'manifest-owned design seed selector missing')
        return fn(dict(construction=dict(design_operation=expression),parts=parts))

    def test_inner_domain_not_outer_surface(self):
        self.assertEqual(self.select('DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',
            {'outer-ellipsoid':{},'inner-ellipsoid-candidate':{}}),'inner-ellipsoid-candidate')

    def test_missing_seed_cannot_fall_back_to_outer(self):
        with self.assertRaises(ValueError):
            self.select('DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',{'outer-ellipsoid':{}})

    def test_unknown_formula_cannot_guess_a_domain(self):
        with self.assertRaises(ValueError):self.select('UNION(outer-ellipsoid, ALL fixed_parts)',{'outer-ellipsoid':{}})
