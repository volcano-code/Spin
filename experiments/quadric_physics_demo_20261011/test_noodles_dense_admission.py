import copy
import json
from pathlib import Path
import unittest


class DenseAdmissionTests(unittest.TestCase):
    def test_versioned_identity_cannot_exchange_inputs(self):
        from noodles_admission import validate_manifest
        here=Path(__file__).resolve().parent
        dense=json.loads((here/'noodles-dense-v2/manifest.json').read_text())
        old=json.loads((here/'noodles-input-v1/manifest.json').read_text())
        validate_manifest(dense,'dense-v2')
        validate_manifest(old,'v1')
        with self.assertRaises(ValueError):validate_manifest(dense,'v1')
        with self.assertRaises(ValueError):validate_manifest(old,'dense-v2')
        changed=copy.deepcopy(dense);changed['construction']['fixed_parts'].pop()
        with self.assertRaisesRegex(ValueError,'FIXED_PARTS'):validate_manifest(changed,'dense-v2')

    def test_overlap_screen_does_not_claim_a_boolean_proof(self):
        from noodles_dense_admission import resource_and_overlap_screen
        here=Path(__file__).resolve().parent
        result=resource_and_overlap_screen(here/'noodles-dense-v2')
        self.assertEqual(result['status'],'ESTIMATE_AND_BROAD_PHASE_ONLY')
        self.assertGreater(result['tube_pair_aabb_overlap_count'],0)
        self.assertEqual(result['planned_exact_boolean_calls'],325)
        self.assertFalse(result['whole_union_validated'])
        self.assertFalse(result['solver_admitted'])


if __name__=='__main__':unittest.main()
