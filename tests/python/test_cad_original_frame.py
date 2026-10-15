"""Restore physical coordinates, not merely a mass number in a sidecar."""
import json
import copy
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from test_exact_material_region import box

ROOT=Path(__file__).resolve().parents[2]
PROBE=Path(os.environ['MAKE_IT_STAND_OCCT_STEP_PROBE'])


class OriginalFrameCadTests(unittest.TestCase):
    def test_actual_shifted_multidensity_b_keeps_physical_mass_and_frame(self):
        from make_it_stand.verification.cad_coordinates import export_original_frame_cad
        from make_it_stand.verification.buoyant_cad import export_saved_buoyant_step
        from make_it_stand.verification.spin_bundle import load_spin_artifact
        from test_buoyant_table_dispatch import buoyant_case
        from test_spin_it_faster_table1_runner import load_runner
        case=buoyant_case()
        p=case['buoyant_problem']; offset=np.array([3.,-4.,5.])
        for key in ('outer','design','displacement'):
            p[key]['vertices']=(2*np.array(p[key]['vertices'])+offset).tolist()
        for region in p['fixed_regions']:
            region['vertices']=(2*np.array(region['vertices'])+offset).tolist()
            region['density']*=3
        p['design']['density']*=3;p['liquid_density']*=3
        p['target_origin']=offset.tolist()
        for stage in (p,p['stage_b']):
            stage['waterline']['origin']=(2*np.array(stage['waterline']['origin'])+offset).tolist()
            stage['initial_parameters']=(np.array(stage['initial_parameters'])/2).tolist()
        p['length_unit']='cm';p['minimum_outer_wall']*=2
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            solved=load_runner().run_case(case,1,1,bundle_root=root/'materials')
            normalized=export_saved_buoyant_step(solved,bundle_root=root/'materials',run_index=0,
                stage='B',occt_executable=PROBE,output_directory=root/'normalized')
            stage=solved['runs'][0]['stages'][1]
            artifact=load_spin_artifact(root/'materials'/stage['relative_bundle'],
                                       expected_problem_digest=stage['problem_digest'])
            report=export_original_frame_cad(artifact,normalized['cad'],
                normalized_brep=root/'normalized/material.brep',occt_executable=PROBE,
                output_directory=root/'original')
            self.assertEqual(report['length_unit'],'cm')
            self.assertEqual(report['material_region_ids'],['design:B','A:shell','A:design:A'])
            self.assertAlmostEqual(report['raw_mass_moments']['step'][0],737.28,places=7)
            np.testing.assert_allclose(np.array(report['raw_mass_moments']['step'][1:3])/737.28,
                                       [3.,-4.],rtol=0,atol=1e-9)
            self.assertLessEqual(report['maximum_normalized_moment_error'],1e-6)
            self.assertTrue((root/'original/material.step').is_file())
            for key in ('entity_densities','entity_solid_counts','artifact_scientific_digest','material_entity_ids'):
                forged=copy.deepcopy(normalized['cad'])
                forged[key]=[9]*len(forged[key]) if isinstance(forged[key],list) else '0'*64
                with self.subTest(field=key),self.assertRaises(ValueError):
                    export_original_frame_cad(artifact,forged,normalized_brep=root/'normalized/material.brep',
                        occt_executable=PROBE,output_directory=root/key)

    def test_malformed_native_response_has_a_failure_record_not_a_published_cad(self):
        from make_it_stand.verification.cad_coordinates import export_original_frame_cad
        from make_it_stand.verification.buoyant_cad import export_saved_buoyant_step
        from make_it_stand.verification.spin_bundle import load_spin_artifact
        from test_buoyant_table_dispatch import buoyant_case
        from test_spin_it_faster_table1_runner import load_runner
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            solved=load_runner().run_case(buoyant_case(),1,1,bundle_root=root/'materials')
            normalized=export_saved_buoyant_step(solved,bundle_root=root/'materials',run_index=0,
                stage='B',occt_executable=PROBE,output_directory=root/'normalized')
            stage=solved['runs'][0]['stages'][1]
            artifact=load_spin_artifact(root/'materials'/stage['relative_bundle'],
                                       expected_problem_digest=stage['problem_digest'])
            response=subprocess.CompletedProcess(['occt'],0,stdout='[]',stderr='')
            with patch('make_it_stand.verification.cad_coordinates.subprocess.run',return_value=response):
                with self.assertRaises(ValueError):
                    export_original_frame_cad(artifact,normalized['cad'],
                        normalized_brep=root/'normalized/material.brep',occt_executable=PROBE,
                        output_directory=root/'failed')
            self.assertEqual(json.loads((root/'failed/failure.json').read_text())['type'],'ValueError')
            self.assertFalse((root/'failed/material.step').exists())
            self.assertFalse((root/'failed/material.step.identity.json').exists())

    def test_scaled_translated_solid_and_step_units_are_real_geometry(self):
        # Unit cube -> side 2, center (3,-4,5), in metres.
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            mesh=box((-.5,)*3,(.5,)*3)
            obj=root/'box.obj'
            obj.write_text(''.join('v %.17g %.17g %.17g\n'%tuple(v) for v in mesh[0])+
                           ''.join('f %d %d %d\n'%tuple(f+1) for f in mesh[1]))
            source=subprocess.run([str(PROBE),'--mesh',str(obj),'1',str(root/'n.brep'),str(root/'n.step')],
                                  capture_output=True,text=True,timeout=60)
            self.assertEqual(source.returncode,0,source.stderr)
            done=subprocess.run([str(PROBE),'--restore-material-frame',str(root/'n.brep'),
                '2','3','-4','5','m',str(root/'p.brep'),str(root/'p.step')],
                capture_output=True,text=True,timeout=60)
            self.assertEqual(done.returncode,0,done.stderr)
            report=json.loads(done.stdout)
            want=np.array([8.,24.,-32.,40.,224/3,392/3,608/3,-160.,120.,-96.])
            for name in ('transformed','brep','step'):
                np.testing.assert_allclose(report['solid_geometric_moments'][name],[want],rtol=1e-10,atol=1e-10)
            # Independently reading STEP in mm catches a wrong/missing file unit.
            np.testing.assert_allclose(report['solid_geometric_moments']['step_mm'],
                [want*np.array([1e9,*([1e12]*3),*([1e15]*6)])],rtol=1e-10,atol=1e-10)
            self.assertEqual(report['length_unit'],'m')
            with self.subTest('no overwrite'):
                before=(root/'p.step').read_bytes()
                repeated=subprocess.run([str(PROBE),'--restore-material-frame',str(root/'n.brep'),
                    '2','3','-4','5','m',str(root/'p.brep'),str(root/'p.step')],capture_output=True,text=True)
                self.assertNotEqual(repeated.returncode,0)
                self.assertEqual((root/'p.step').read_bytes(),before)


if __name__=='__main__':unittest.main()
