"""Real fixed-normal source and reference Boolean; no optimizer or paper grant."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from make_it_stand import spin
import test_faster_cherchi_basis as shared

FIXTURE=shared.ROOT/'experiments/spin_it_reproduction/faster/fixtures/cherchi-fixed-normal-paraboloid-method-v1/problem.json'
MAP='BOUNDED_FINITE_NORMAL_PARABOLOID_V1'

class FixedNormalParaboloidTests(unittest.TestCase):
    def test_native_interface_retains_explicit_source_and_real_eight_dimensions(self):
        native=spin._native_module()
        settings=native.SpinTriangulationSettings()
        self.assertTrue(hasattr(settings,'fixed_normal_paraboloid'),'fixed-normal Python request missing')
        self.assertFalse(settings.fixed_normal_paraboloid)
        p=json.loads(FIXTURE.read_text())
        settings.icosphere_subdivisions=2;settings.truncation_radius=16.
        settings.paraboloid_grid_resolution=16;settings.maximum_radius_doublings=0
        legacy=native.evaluate_spin_it_faster_stratum_source_motion(p['coefficients'],p['design_vertices'],settings)
        settings.fixed_normal_paraboloid=True
        source=native.evaluate_spin_it_faster_stratum_source_motion(p['coefficients'],p['design_vertices'],settings)
        self.assertTrue(source['defined'],source['diagnostic'])
        self.assertTrue(source['fixed_normal_paraboloid'])
        self.assertEqual(np.asarray(source['coefficient_basis']).shape,(10,8))
        self.assertEqual(np.asarray(source['vertex_velocity']).shape,(438,8))
        self.assertEqual(np.asarray(source['fixed_normals']).shape,(162,3))
        self.assertEqual(len(source['source_normal_vertices']),146)
        self.assertNotEqual(source['triangulation_hash'],legacy['triangulation_hash'])
        full=native.evaluate_spin_it_faster_source_motion(p['coefficients'],p['design_vertices'],settings)
        self.assertFalse(full['defined'])

    def test_invalid_profile_cannot_be_silently_ignored(self):
        runner=shared.CherchiBasisTests().module()
        native=spin._native_module();p=json.loads(FIXTURE.read_text())
        for value in (1,'true',None):
            with self.subTest(value=value),self.assertRaises(ValueError):
                runner._settings(native,dict(p,fixed_normal_paraboloid=value))
        for patch in ({'paper_projective_hyperboloid':True},{'algebraic_finite_map':True},
                      {'schema':'cherchi-full-basis-method-v1'},{'maximum_icosphere_subdivisions':3}):
            with self.subTest(patch=patch),self.assertRaises(ValueError):
                runner._settings(native,dict(p,**patch))

    def test_actual_fixed_normal_material_checks_all_eight_directions(self):
        runner=shared.CherchiBasisTests().module();p=json.loads(FIXTURE.read_text())
        self.assertEqual(runner._protocol(p).get('source_map'),MAP,'request selected the legacy grid')
        with tempfile.TemporaryDirectory() as temporary:
            out=Path(temporary)/'out'
            report=runner.run(out,fixture=FIXTURE,probe=shared.PROBE)
            self.assertEqual(report['status'],'PASSED',json.dumps(report))
            self.assertEqual(report['evaluation_count'],49)
            self.assertEqual(len(report['checks']),8)
            self.assertFalse(report['optimizer_executed'])
            rows=[json.loads(s) for s in (out/'evaluations.jsonl').read_text().splitlines()]
            self.assertEqual(np.asarray(rows[0]['jacobian']).shape,(10,8))
            self.assertTrue(all(r['independent_moment_pass'] for r in rows))
            self.assertIn('LPI',rows[0]['provenance_types'])
            for defect in ('missing','changed','legacy','padded','surface_type'):
                bad=copy.deepcopy(rows)
                if defect=='missing':bad[0].pop('fixed_normal_source_features')
                elif defect=='changed':bad[0]['fixed_normal_source_features']['normal_vertices'][0]+=1
                elif defect=='legacy':bad[0]['branch']['source_map']='LEGACY_GRID'
                elif defect=='padded':bad[0]['jacobian']=np.pad(bad[0]['jacobian'],((0,0),(0,2))).tolist()
                else:
                    for row in bad: row['branch']['surface_type']='HYPERBOLOID'
                self.assertNotEqual(runner.capture_comparison(report['protocol'],bad)['status'],'PASSED')
            for source_map in (None, 'LEGACY_GRID', 'UNKNOWN_SOURCE'):
                downgraded=copy.deepcopy(report['protocol'])
                if source_map is None: downgraded.pop('source_map')
                else: downgraded['source_map']=source_map
                self.assertNotEqual(runner.capture_comparison(downgraded,rows)['status'],'PASSED',
                                    'fixed-normal records accepted under a downgraded protocol')
            empty=copy.deepcopy(report['protocol']);empty['directions']=[]
            self.assertNotEqual(runner.capture_comparison(empty,rows)['status'],'PASSED',
                                'no checked directions cannot pass a fixed-normal stencil')

if __name__=='__main__':unittest.main()
