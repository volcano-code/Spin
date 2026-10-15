"""Actual saved-bundle output path, never an optimization or manufacture test."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import test_derived_spin_material as fixtures
from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact

ROOT=Path(__file__).resolve().parents[2]
RUNNER=ROOT/'experiments/spin_it_reproduction/strict/derive_saved_spin_material.py'


class DerivedMaterialRunnerTests(unittest.TestCase):
    def test_success_syncs_parent_entry_as_well_as_material_and_receipt(self):
        import importlib.util
        from unittest.mock import patch
        spec=importlib.util.spec_from_file_location('derived_output_runner',RUNNER)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        authority=fixtures.DerivedMaterialTests().fixture()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'authority';save_spin_artifact(authority,source)
            parent=root.stat();synced=[];real_fsync=os.fsync
            def observe_and_sync(fd):
                item=os.fstat(fd)
                real_fsync(fd)  # Observe actual fsync, not a mocked write success.
                synced.append((item.st_dev,item.st_ino))
            with patch('os.fsync',side_effect=observe_and_sync):
                result=module.run(source,root/'derived',
                    expected_problem_digest=authority.problem.input_digest,
                    expected_scientific_digest=authority.scientific_digest)
            self.assertEqual(result['relation']['derived_scientific_digest'],authority.scientific_digest)
            self.assertIn((parent.st_dev,parent.st_ino),synced,
                          'successful export did not fsync its parent directory entry')

    def invoke(self,source,output,authority,**kw):
        self.assertTrue(RUNNER.is_file(),'formal saved-material output runner is missing')
        return subprocess.run([sys.executable,str(RUNNER),'--source',str(source),
            '--output',str(output),'--expected-problem-digest',authority.problem.input_digest,
            '--expected-scientific-digest',kw.get('digest',authority.scientific_digest)],
            capture_output=True,text=True,timeout=30,env=os.environ.copy())

    def test_independent_processes_recheck_persisted_geometry_and_preserve_history(self):
        authority=fixtures.DerivedMaterialTests().fixture()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'authority';save_spin_artifact(authority,source)
            before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir()}
            receipts=[]
            for name in ('first','repeat'):
                output=root/name;run=self.invoke(source,output,authority)
                self.assertEqual(run.returncode,0,run.stderr)
                record=json.loads((output/'derivation.json').read_bytes())
                self.assertEqual(record['source_reaudit']['source_manifest_sha256'],before['manifest.json'])
                self.assertEqual(record['source_reaudit']['solver_freshness'],'NOT_REASSESSED')
                self.assertEqual(record['relation']['authority_scientific_digest'],authority.scientific_digest)
                loaded=load_spin_artifact(output/'material',expected_problem_digest=authority.problem.input_digest)
                self.assertEqual(record['relation']['derived_scientific_digest'],loaded.scientific_digest)
                self.assertEqual(loaded.facts['physical']['mass'],7/8)
                self.assertEqual(loaded.facts['manufacturing_admission'],'BLOCKED_UNKNOWN_WALL_THICKNESS')
                self.assertEqual(record['material_manifest']['sha256'],hashlib.sha256(
                    (output/'material/manifest.json').read_bytes()).hexdigest())
                receipts.append(record)
            self.assertEqual(receipts[0],receipts[1])
            self.assertEqual(before,{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir()})
            duplicate=self.invoke(source,root/'first',authority)
            self.assertNotEqual(duplicate.returncode,0)
            self.assertEqual(json.loads((root/'first/derivation.json').read_bytes()),receipts[0])

    def test_wrong_source_identity_and_corrupted_saved_payload_cannot_publish_success(self):
        authority=fixtures.DerivedMaterialTests().fixture()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'authority';save_spin_artifact(authority,source)
            run=self.invoke(source,root/'wrong',authority,digest='0'*64)
            self.assertNotEqual(run.returncode,0)
            self.assertFalse((root/'wrong/derivation.json').exists())
            mesh=next(source.glob('*vertices.npy'))
            mesh.write_bytes(mesh.read_bytes()+b'altered')
            run=self.invoke(source,root/'corrupt',authority)
            self.assertNotEqual(run.returncode,0)
            self.assertFalse((root/'corrupt/derivation.json').exists())


if __name__=='__main__':unittest.main()
