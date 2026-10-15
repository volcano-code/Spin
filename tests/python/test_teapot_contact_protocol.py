"""Synthetic protocol regressions; these are not physical experiment evidence."""
import copy
import importlib.util
import math
from pathlib import Path
import sys
import unittest

path = Path(__file__).resolve().parents[2] / 'experiments/teapot_spin2014/contact_protocol.py'
spec = importlib.util.spec_from_file_location('teapot_contact_protocol', path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def fixture():
    rows = []
    for repeat in [0, 1]:
        for profile in ['default', 'robust']:
            for rate in [240,480,960]:
                for d in [None,0,1,2,3,4,5,6,7]:
                    tilt = 0.0 if d is None else 2.0
                    rows.append(dict(repeat=repeat,profile=profile,step_denominator=rate,direction=d,
                        initial_tilt=tilt,azimuth=0 if d is None else d*math.pi/4,
                        dt=1/rate,duration=8.,hz=35.,friction=.8,artifact_digest='a'*64,
                        protocol='PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT',constraint_count=0,
                        samples=[dict(time=0.,tilt_degrees=tilt),dict(time=8.,tilt_degrees=tilt)]))
    return rows


class ProtocolTests(unittest.TestCase):
    def audit(self, rows):
        return module.audit_records(rows,artifact_digest='a'*64,spin_hz=35.)

    def test_complete_coverage_does_not_grant_physics_or_determinism(self):
        result = self.audit(fixture())
        self.assertEqual((result['observed_scenarios'],result['observed_nominal']),(108,12))
        self.assertEqual(result['status'],'COMPLETE_PROTOCOL_COVERAGE')
        self.assertEqual(result['physics_acceptance'],'NOT_GRANTED')
        self.assertEqual(result['determinism'],'NOT_EVALUATED')

    def test_direction_zero_cannot_be_called_nominal(self):
        rows=fixture();rows[0]['initial_tilt']=2.;rows[0]['samples'][0]['tilt_degrees']=2.
        with self.assertRaisesRegex(ValueError,'initial_tilt'): self.audit(rows)

    def test_missing_repeat_and_nominal_cannot_be_hidden(self):
        rows=[r for r in fixture() if r['repeat']==0 and r['direction'] is not None]
        result=self.audit(rows)
        self.assertEqual(result['status'],'INCOMPLETE')
        self.assertEqual(result['observed_nominal'],0)
        self.assertEqual(len(result['missing_scenarios']),60)

    def test_duplicate_wrong_artifact_nan_and_step_mismatch_rejected(self):
        rows=fixture()
        with self.assertRaisesRegex(ValueError,'DUPLICATE'):self.audit(rows+[copy.deepcopy(rows[0])])
        for key,value,reason in [('artifact_digest','b'*64,'ARTIFACT'),('dt',1/2000,'dt'),
                                 ('hz',float('nan'),'NONFINITE'),('constraint_count',1,'FREE_CONTACT')]:
            bad=copy.deepcopy(rows);bad[0][key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,reason):self.audit(bad)

    def test_fake_pass_and_truncated_trajectory_do_not_pass(self):
        rows=fixture()[:1];rows[0]['passed']=True
        self.assertEqual(self.audit(rows)['status'],'INCOMPLETE')
        rows[0]['samples'][-1]['time']=7.
        with self.assertRaisesRegex(ValueError,'sample_end'):self.audit(rows)

    def test_missing_or_boolean_identity_cannot_be_nominal(self):
        for field, value in [('direction',False),('repeat',False),('profile',{}),('step_denominator',True)]:
            row=fixture()[0];row[field]=value
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'EXPLICIT_SCENARIO'):
                self.audit([row])
        row=fixture()[0];del row['direction']
        with self.assertRaisesRegex(ValueError,'EXPLICIT_SCENARIO'):self.audit([row])


if __name__=='__main__':unittest.main()
