"""Actual Cherchi plane map: three directions, never a padded top Jacobian."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import test_faster_cherchi_basis as existing

FIXTURE = existing.ROOT/'experiments/spin_it_reproduction/faster/fixtures/cherchi-plane-basis-method-v1/problem.json'


class CherchiPlaneBasisTests(unittest.TestCase):
    def test_plane_protocol_preserves_only_three_tangent_directions(self):
        runner = existing.CherchiBasisTests().module()
        problem = json.loads(FIXTURE.read_text())
        protocol = runner._protocol(problem)
        basis = np.asarray(protocol['parameter_basis'])
        q = np.asarray(protocol['normalized_coefficients'])
        self.assertEqual(basis.shape, (10, 3))
        np.testing.assert_allclose(basis.T@basis, np.eye(3), atol=1e-13)
        np.testing.assert_allclose(q@basis, 0, atol=1e-13)
        np.testing.assert_array_equal(basis[4:], 0)
        self.assertEqual(len(protocol['samples']), 19)
        for row in protocol['samples']:
            np.testing.assert_array_equal(np.asarray(row['coefficients'])[4:], 0)
        curved = copy.deepcopy(problem); curved['coefficients'][8] = 1e-8
        with self.assertRaisesRegex(ValueError, 'plane|Plane'):
            runner._protocol(curved)

    def test_actual_plane_boolean_all_directions_and_independent_moments(self):
        runner = existing.CherchiBasisTests().module()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'output'
            report = runner.run(output, fixture=FIXTURE, probe=existing.PROBE)
            self.assertEqual(report['status'], 'PASSED', json.dumps(report))
            self.assertFalse(report['optimizer_executed'])
            self.assertEqual(report['evaluation_count'], 19)
            self.assertEqual(len(report['checks']), 3)
            self.assertIn('10x3', report['scope'])
            rows = [json.loads(x) for x in (output/'evaluations.jsonl').read_text().splitlines()]
            self.assertEqual(np.asarray(rows[0]['jacobian']).shape, (10,3))
            self.assertTrue(all(r['independent_moment_pass'] for r in rows))
            self.assertTrue(all(r['branch']['surface_type'] == 'PLANE' for r in rows))
            self.assertIn('LPI', rows[0]['provenance_types'])
            for defect in ('padded', 'column', 'missing', 'branch', 'chart'):
                bad = copy.deepcopy(rows)
                if defect == 'padded':
                    bad[0]['jacobian'] = np.pad(bad[0]['jacobian'], ((0,0),(0,7))).tolist()
                elif defect == 'column': bad[0]['jacobian'][0][2] += .5
                elif defect == 'missing': bad.pop()
                elif defect == 'branch': bad[1]['branch']['surface_type'] = 'PARABOLOID'
                else: bad[1]['branch']['plane_frame_chart'] = 1
                with self.subTest(defect=defect):
                    self.assertNotEqual(runner.capture_comparison(report['protocol'], bad)['status'], 'PASSED')
            for defect in ('missing_direction', 'empty_directions', 'empty_steps', 'wrong_basis', 'changed_sample'):
                protocol = copy.deepcopy(report['protocol'])
                if defect == 'missing_direction': protocol['directions'].pop()
                elif defect == 'empty_directions': protocol['directions'] = []
                elif defect == 'empty_steps': protocol['steps'] = []
                elif defect == 'wrong_basis': protocol['parameter_basis'][8][0] = .1
                else:
                    # Even changing both copies must not bypass the frozen stencil.
                    bad = copy.deepcopy(rows)
                    protocol['samples'][1]['coefficients'][0] += .01
                    bad[1]['coefficients'] = protocol['samples'][1]['coefficients']
                    self.assertNotEqual(runner.capture_comparison(protocol,bad)['status'], 'PASSED')
                    continue
                with self.subTest(protocol_defect=defect):
                    self.assertNotEqual(runner.capture_comparison(protocol,rows)['status'], 'PASSED')
            absent_chart = copy.deepcopy(rows)
            for row in absent_chart: row['branch'].pop('plane_frame_chart',None)
            self.assertNotEqual(runner.capture_comparison(report['protocol'],absent_chart)['status'], 'PASSED')

    def test_native_plane_frame_seam_is_explicit_not_a_reused_derivative(self):
        runner = existing.CherchiBasisTests().module()
        problem = json.loads(FIXTURE.read_text())
        for nx, chart in ((.9-1e-5,0),(.9+1e-5,1)):
            problem['coefficients'] = [.17,nx,float(np.sqrt(1-nx*nx)),0,0,0,0,0,0,0]
            protocol = runner._protocol(problem)
            self.assertEqual(protocol['plane_frame_chart'],chart)
        problem['coefficients'] = [.17,.9,float(np.sqrt(1-.9**2)),0,0,0,0,0,0,0]
        with self.assertRaisesRegex(ValueError,'chart|boundary'):
            runner._protocol(problem)

    def test_actual_chart_crossing_does_not_pass_fixed_branch_audit(self):
        runner = existing.CherchiBasisTests().module()
        problem = json.loads(FIXTURE.read_text())
        nx = .9-1e-5
        problem['coefficients'] = [.17,nx,float(np.sqrt(1-nx*nx)),0,0,0,0,0,0,0]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            fixture = directory/'seam-neighborhood.json'
            fixture.write_text(json.dumps(problem))
            report = runner.run(directory/'output',fixture=fixture,probe=existing.PROBE)
            rows = [json.loads(x) for x in (directory/'output/evaluations.jsonl').read_text().splitlines()]
            self.assertTrue(all(r['status'] == 'VALID_EXACT_MATERIAL' for r in rows))
            self.assertEqual({r['branch']['plane_frame_chart'] for r in rows},{0,1})
            self.assertEqual(report['status'],'FAILED_DIFFERENTIAL_OR_BRANCH')
            self.assertTrue(any(not s['same_branch'] for c in report['checks'] for s in c['steps']))


if __name__ == '__main__':
    unittest.main()
