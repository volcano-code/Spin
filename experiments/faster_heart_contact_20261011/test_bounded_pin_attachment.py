import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


class AttachmentTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('bounded_pin_attachment'))
        import bounded_pin_attachment
        return bounded_pin_attachment

    def test_unknown_or_source_mismatch_cannot_be_upgraded_by_true_field(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'terminal.json'
            path.write_text(json.dumps(dict(schema='bounded-pin-companion-admission-v1',
                status='PIN_CORRESPONDENCE_UNKNOWN',pin_correspondence_certified=True)))
            with self.assertRaisesRegex(ValueError,'NOT_CERTIFIED'):
                self.module().load_companion(path,expected_source_hashes={})

    def test_qualified_label_without_actual_files_and_bound_sources_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'terminal.json'
            path.write_text(json.dumps(dict(schema='bounded-pin-companion-admission-v1',
                status='PIN_CORRESPONDENCE_CERTIFIED',pin_correspondence_certified=True,
                inputs={},files={})))
            with self.assertRaisesRegex(ValueError,'SOURCE_OR_FILES_REQUIRED'):
                self.module().load_companion(path,expected_source_hashes={})

    def test_recomputed_identity_correspondence_does_not_grant_whole_request(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);out=root/'out';out.mkdir()
            pin=root/'pin.obj';pin.write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\n'
                'f 1 3 2\nf 1 2 4\nf 1 4 3\nf 2 3 4\n')
            mesh=dict(vertices=[['0/1']*3,['1/1','0/1','0/1'],['0/1','1/1','0/1'],['0/1','0/1','1/1']],
                      faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
            def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
            outer=root/'outer.json';outer.write_text(json.dumps(mesh))
            fixed=root/'fixed.json';fixed.write_text(json.dumps(mesh))
            inputs={name:dict(path=str(path),sha256=digest(path)) for name,path in
                [('pin',pin),('authority-outer',outer),('candidate-fixed.rational.json',fixed)]}
            derived=out/'derived-pin.rational.json'
            derived.write_text(json.dumps(dict(mesh,role='BOUNDED_DERIVED_PIN_NOT_ORIGINAL_PIN',
                construction='ORIGINAL_PIN_INTERSECT_SAME_DERIVED_FIXED',original_pin_sha256=digest(pin),
                derived_fixed_sha256=digest(fixed))))
            identity=out/'identity.json';identity.write_text(json.dumps(dict(normalization_length_authority_exact='1/1')))
            terminal=out/'terminal.json';terminal.write_text(json.dumps(dict(schema='bounded-pin-companion-admission-v1',
                status='PIN_CORRESPONDENCE_CERTIFIED',pin_correspondence_certified=True,
                actual_removed_volume_exact='0/1',inputs=inputs,
                files={p.name:digest(p) for p in (derived,identity)})))
            report=self.module().load_companion(terminal,expected_source_hashes={k:v['sha256'] for k,v in inputs.items()})
            self.assertTrue(report['native_fresh_checks_required'])
            self.assertFalse(report['geometry_admission_granted'])
            derived.write_text(derived.read_text()+' ')
            with self.assertRaisesRegex(ValueError,'FILE_IDENTITY'):
                self.module().load_companion(terminal,expected_source_hashes={k:v['sha256'] for k,v in inputs.items()})

if __name__=='__main__':unittest.main()
