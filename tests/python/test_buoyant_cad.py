"""Actual A→B material must reach CAD without changing its run identity."""
import copy
import json
import os
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from test_buoyant_table_dispatch import buoyant_case
from test_spin_it_faster_table1_runner import load_runner
from make_it_stand.verification.spin_common import digest

ROOT=Path(__file__).resolve().parents[2]
PROBE=Path(os.environ.get('MAKE_IT_STAND_OCCT_STEP_PROBE',
    ROOT/'build/occt-buoyant-chain-20261001/occt_step_roundtrip'))


class BuoyantCadTests(unittest.TestCase):
    def test_outer_composition_cannot_contradict_saved_a_and_b(self):
        from make_it_stand.verification.buoyant_cad import export_saved_buoyant_step
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result=load_runner().run_case(buoyant_case(),1,1,bundle_root=root/'materials')
            changes=[(key,'0'*64) for key in ('source_a_scientific_digest',
                'source_a_evidence_digest','artifact_scientific_digest',
                'b_artifact_scientific_digest','b_artifact_evidence_digest')]
            changes += [('status','FAILED'),('solver_status','FAILED'),('seal_error','not sealed')]
            for key,value in changes:
                with self.subTest(field=key):
                    forged=copy.deepcopy(result)
                    forged['runs'][0]['stages'][1]['seal_identity'][key]=value
                    forged['scientific_hash']=digest({k:v for k,v in forged.items() if k!='scientific_hash'})
                    with patch('make_it_stand.verification.buoyant_cad.export_spin_it_faster_quadric_material_step',
                               side_effect=AssertionError('inconsistent identity reached CAD')):
                        with self.assertRaises(ValueError):
                            export_saved_buoyant_step(forged,bundle_root=root/'materials',run_index=0,
                                stage='B',occt_executable=PROBE,output_directory=root/key)
                    self.assertFalse((root/key).exists())

    def test_cad_timeout_preserves_original_failure(self):
        from make_it_stand.verification.buoyant_cad import export_saved_buoyant_step
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result=load_runner().run_case(buoyant_case(),1,1,bundle_root=root/'materials')
            error=subprocess.TimeoutExpired([str(PROBE)],60)
            with patch('make_it_stand.verification.faster_material.subprocess.run',side_effect=error):
                with self.assertRaises(subprocess.TimeoutExpired) as caught:
                    export_saved_buoyant_step(result,bundle_root=root/'materials',run_index=0,
                        stage='B',occt_executable=PROBE,output_directory=root/'timeout')
            self.assertIs(caught.exception,error)
            self.assertEqual(json.loads((root/'timeout/failure.json').read_text()),
                             {'type':'TimeoutExpired','message':str(error)})
            self.assertTrue((root/'timeout/request.json').is_file())
            self.assertFalse((root/'timeout/buoyant-cad.json').exists())

    def test_actual_b_exports_its_shell_and_retained_a_without_relabelling(self):
        from make_it_stand.verification.buoyant_cad import export_saved_buoyant_step
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result=load_runner().run_case(buoyant_case(),1,1,bundle_root=root/'materials')
            before=copy.deepcopy(result)
            report=export_saved_buoyant_step(result,bundle_root=root/'materials',run_index=0,
                stage='B',occt_executable=PROBE,output_directory=root/'cad')
            self.assertEqual(result,before)
            self.assertEqual(report['source_operation'],'FLOATING_ISLAND_BUOYANT_B_FROM_SEALED_A')
            self.assertEqual(report['source_a_scientific_digest'],result['runs'][0]['stages'][0]['artifact_scientific_digest'])
            self.assertEqual(report['artifact_scientific_digest'],result['runs'][0]['stages'][1]['artifact_scientific_digest'])
            self.assertEqual(report['material_region_ids'],['design:B','A:shell','A:design:A'])
            self.assertLessEqual(report['maximum_normalized_material_moment_error'],1e-6)
            self.assertAlmostEqual(report['physical_mass'],30.72,places=8)
            for name in ('material.brep','material.step','material.step.identity.json','buoyant-cad.json'):
                self.assertTrue((root/'cad'/name).is_file(),name)
            with self.assertRaises(FileExistsError):
                export_saved_buoyant_step(result,bundle_root=root/'materials',run_index=0,
                    stage='B',occt_executable=PROBE,output_directory=root/'cad')
            forged=copy.deepcopy(result)
            forged['runs'][0]['stages'][1]['source_a_scientific_digest']='0'*64
            forged['scientific_hash']=digest({k:v for k,v in forged.items() if k!='scientific_hash'})
            with self.assertRaises(ValueError):
                export_saved_buoyant_step(forged,bundle_root=root/'materials',run_index=0,
                    stage='B',occt_executable=PROBE,output_directory=root/'forged')


if __name__=='__main__':unittest.main()
