"""Report binding tests; synthetic fixtures are not physical certifications."""
import json
from pathlib import Path
import tempfile
import unittest

import flower_result_figure as api
from test_flower_physics_admission import FlowerPhysicsAdmissionTests, digest, write


class FigureStatusTests(unittest.TestCase):
    def fixture(self, root):
        folder,material,legacy,saved=FlowerPhysicsAdmissionTests().setup_saved_records(root)
        out=root/'physics';out.mkdir()
        rows=[]
        for model in ('ideal','flexible'):
            for variant in ('original','optimized'):
                for phase in range(8):
                    for dt in (.0005,.00025):
                        identity=f'{model}-{variant}-phase{phase}-dt{dt}'
                        trajectory=out/(identity+'.npz')
                        trajectory.write_bytes(b'identity-only test fixture, not physical trajectory')
                        rows.append(dict(id=identity,model=model,variant=variant,phase_index=phase,
                            dt_seconds=dt,actual_steps=round(8/dt),trajectory_sha256=digest(trajectory)))
        proof=json.loads(saved.read_text())
        physics=out/'summary.json'
        write(physics,dict(schema='admitted-flower-actual-physics-matrix-v1',models=['ideal','flexible'],
            phase_count=8,time_steps_seconds=[.0005,.00025],duration_seconds=8.,
            expected_scenes=64,recorded_scenes=64,scenes=rows,motor_drive=False,
            original=proof['original'],optimized=proof['derived_optimized'],
            sources=dict(material_preparation_sha256=digest(material/'output-preparation.json'),
                saved_material_check_path=str(saved),saved_material_check_sha256=digest(saved))))
        return folder,saved,physics

    def validate(self,*args,**kwargs):
        self.assertTrue(callable(getattr(api,'validate_status_sources',None)),
                        'figure must verify stored material and physical records, not hand-label status')
        return api.validate_status_sources(*args,**kwargs)

    def test_complete_bound_sources_enable_physics_status(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,saved,physics=self.fixture(Path(directory))
            result=self.validate(folder,saved_material_check=saved,physics_summary=physics)
            self.assertTrue(result['physics_result'])
            self.assertEqual(result['recorded_scenes'],64)
            self.assertEqual(result['saved_material_check_sha256'],digest(saved))

    def test_mismatched_physics_material_source_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,saved,physics=self.fixture(Path(directory))
            report=json.loads(physics.read_text());report['sources']['material_preparation_sha256']='wrong'
            write(physics,report)
            with self.assertRaisesRegex(ValueError,'PHYSICS_MATERIAL_IDENTITY'):
                self.validate(folder,saved_material_check=saved,physics_summary=physics)

    def test_incomplete_duplicate_or_modified_trajectory_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,saved,physics=self.fixture(Path(directory))
            original=json.loads(physics.read_text())
            for rows in (original['scenes'][:-1],original['scenes'][:-1]+[original['scenes'][0]]):
                write(physics,{**original,'scenes':rows})
                with self.assertRaisesRegex(ValueError,'PHYSICS_SCENE'):
                    self.validate(folder,saved_material_check=saved,physics_summary=physics)
            write(physics,original)
            (physics.parent/(original['scenes'][0]['id']+'.npz')).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'PHYSICS_TRAJECTORY'):
                self.validate(folder,saved_material_check=saved,physics_summary=physics)

    def test_physics_without_material_gate_rejected_and_legacy_stays_unapproved(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,saved,physics=self.fixture(Path(directory))
            with self.assertRaisesRegex(ValueError,'SAVED_MATERIAL_CHECK_REQUIRED'):
                self.validate(folder,physics_summary=physics)
            self.assertFalse(self.validate(folder)['physics_result'])

    def test_stale_saved_application_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,saved,physics=self.fixture(Path(directory))
            application=folder/'solve-job-01/output/application.json'
            application.write_text(application.read_text()+' ')
            with self.assertRaisesRegex(ValueError,'SOURCE_IDENTITY'):
                self.validate(folder,saved_material_check=saved,physics_summary=physics)

if __name__=='__main__':unittest.main()
