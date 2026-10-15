import copy
import importlib.util
from pathlib import Path
import unittest

from make_it_stand.verification.quadric_point_jobs import freeze_point_manifest


class SavedProductionReviewTests(unittest.TestCase):
    def check(self, manifest, records, jacobian):
        path = Path(__file__).with_name('review_saved_production.py')
        self.assertTrue(path.exists(), 'frozen production review is not implemented')
        spec = importlib.util.spec_from_file_location('review', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.checked_rows(manifest, records, jacobian)

    def fixture(self):
        manifest = freeze_point_manifest([1.]*10, [.01]*10,
            sources={k:'a'*64 for k in ('input','native','runner','checker')},
            mapping_settings=dict(coefficient_policy='COEFFICIENT_MPFR1024_V1', precision_bits=1024, normal_level=2))
        records = [dict(index=i, status='INDEPENDENTLY_REINTEGRATED',
            geometry_authority=['1']*10, moments=['1']*10,
            branch=dict(provenance_topology_hash='base',surface_type=1),
            retriangulation=dict(proved=True)) for i in range(81)]
        return manifest, records, [['0']*10 for _ in range(10)]

    def test_exact_witness_only_allows_the_proven_retriangulation_change(self):
        m,r,j = self.fixture()
        r[1]['branch']['provenance_topology_hash'] = 'different-triangulation'
        rows = self.check(m,r,j)
        self.assertEqual(len(rows),30)
        self.assertTrue(all(x['passed'] for x in rows))
        self.assertEqual(rows[0]['status'],'EXACT_RETRIANGULATION_WITNESS')
        r[1]['branch']['surface_type'] = 2
        self.assertFalse(self.check(m,r,j)[0]['passed'])

    def test_failed_witness_or_integration_cannot_be_overridden(self):
        m,r,j = self.fixture()
        r[1]['branch']['provenance_topology_hash'] = 'changed'
        r[1]['retriangulation']['proved'] = False
        self.assertFalse(self.check(m,r,j)[0]['passed'])
        r[1]['retriangulation']['proved'] = True
        r[1]['status'] = 'FAILED'
        r[1]['passed'] = True
        self.assertFalse(self.check(m,r,j)[0]['passed'])

    def test_missing_duplicate_points_and_jacobian_error_rejected(self):
        m,r,j = self.fixture()
        with self.assertRaises(ValueError): self.check(m,r[:-1],j)
        duplicate = copy.deepcopy(r); duplicate[1]['index']=0
        with self.assertRaises(ValueError): self.check(m,duplicate,j)
        r[1]['moments'][0]='2'
        self.assertFalse(self.check(m,r,j)[0]['passed'])


if __name__ == '__main__':
    unittest.main()
