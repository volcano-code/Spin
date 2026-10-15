"""Fast admission contract tests, without loading a native module or Boolean."""
import unittest
import gmpy2
import numpy as np
import trimesh
from noodles_admission import boundary_components, translated, measured_interface_clearance


class AdmissionContractTests(unittest.TestCase):
    def test_disconnected_outer_boundary_is_detected(self):
        a = trimesh.creation.box()
        self.assertEqual(boundary_components(a.faces.tolist()), 1)
        faces = a.faces.tolist()+(a.faces+len(a.vertices)).tolist()
        self.assertEqual(boundary_components(faces), 2)

    def test_target_translation_retains_exact_coordinate_information(self):
        region = ([['1/3', '2/3', '-1/6']], [[0,0,0]])
        transformed = translated(region, [0.,0.,-.16])
        self.assertEqual(gmpy2.mpq(transformed[0][0][0]), gmpy2.mpq(1,3))
        self.assertEqual(gmpy2.mpq(transformed[0][0][2]), gmpy2.mpq(-1,6)-gmpy2.mpq(-.16))
        self.assertEqual(transformed[1], region[1])

    def test_observed_isolation_bound_is_conservative_and_rejects_contact(self):
        record = dict(completed=True, applicable=True, squared_distance_exact='1/9')
        bound = measured_interface_clearance(record)
        self.assertGreater(bound,0)
        self.assertLessEqual(gmpy2.mpq(bound)**2, gmpy2.mpq(1,36))
        for changed in (dict(record,squared_distance_exact='0'),dict(record,completed=False),dict(record,applicable=False)):
            with self.assertRaises(ValueError): measured_interface_clearance(changed)


if __name__ == '__main__': unittest.main()
