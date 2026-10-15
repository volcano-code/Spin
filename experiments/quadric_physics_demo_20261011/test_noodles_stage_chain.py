import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


class StageChainTests(unittest.TestCase):
    def make_chain(self, root):
        names = ['hub', 'pin', 'tube-0']
        paths = []
        previous = None
        for index, count in enumerate((2, 3)):
            path = root / str(index)
            path.mkdir()
            boundary = path / 'design.rational.json'
            boundary.write_text(json.dumps({'vertices': [], 'faces': []}))
            packet = dict(schema='noodles-exact-domain-preprocessing-stage-v1',
                manifest_sha256='input', native_sha256='native',
                excluded_parts=names[:count], part_count=3,
                design_sha256=hashlib.sha256(boundary.read_bytes()).hexdigest(),
                predecessor_stage_sha256=previous,
                preprocessing_complete=count == 3,
                region_facts={'completed': True, 'valid': True},
                geometry_admission=False, optimization_admission=False,
                simulation_admission=False)
            record = path / 'stage.json'
            record.write_text(json.dumps(packet))
            previous = hashlib.sha256(record.read_bytes()).hexdigest()
            paths.append(path)
        return names, paths

    def test_complete_chain_only_grants_preprocessing(self):
        import noodles_stage_chain as module
        with TemporaryDirectory() as directory:
            names, paths = self.make_chain(Path(directory))
            report = module.verify_chain(paths, names, 'input', 'native')
        self.assertTrue(report['preprocessing_complete'])
        self.assertFalse(report['geometry_admission'])
        self.assertFalse(report['optimization_admission'])

    def test_missing_stage_and_modified_boundary_rejected(self):
        import noodles_stage_chain as module
        with TemporaryDirectory() as directory:
            names, paths = self.make_chain(Path(directory))
            with self.assertRaisesRegex(ValueError, 'PREDECESSOR'):
                module.verify_chain(paths[1:], names, 'input', 'native')
            (paths[0] / 'design.rational.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'BOUNDARY'):
                module.verify_chain(paths, names, 'input', 'native')

    def test_incomplete_prefix_rejected_for_final_admission(self):
        import noodles_stage_chain as module
        with TemporaryDirectory() as directory:
            names, paths = self.make_chain(Path(directory))
            with self.assertRaisesRegex(ValueError, 'INCOMPLETE'):
                module.verify_chain(paths[:1], names, 'input', 'native')
            report = module.verify_chain(paths[:1], names, 'input', 'native', require_complete=False)
            self.assertEqual(report['excluded_part_count'], 2)

    def test_reordered_exclusions_and_wrong_native_rejected(self):
        import noodles_stage_chain as module
        with TemporaryDirectory() as directory:
            names, paths = self.make_chain(Path(directory))
            with self.assertRaisesRegex(ValueError, 'NATIVE'):
                module.verify_chain(paths, names, 'input', 'other')
            record = paths[0] / 'stage.json'
            packet = json.loads(record.read_text())
            packet['excluded_parts'] = ['pin', 'hub']
            record.write_text(json.dumps(packet))
            with self.assertRaisesRegex(ValueError, 'PREFIX'):
                module.verify_chain(paths, names, 'input', 'native')

    def test_false_completion_and_admission_claims_rejected(self):
        import noodles_stage_chain as module
        with TemporaryDirectory() as directory:
            names, paths = self.make_chain(Path(directory))
            record = paths[0] / 'stage.json'
            packet = json.loads(record.read_text())
            packet['preprocessing_complete'] = True
            record.write_text(json.dumps(packet))
            with self.assertRaisesRegex(ValueError, 'FALSE_COMPLETION'):
                module.verify_chain(paths, names, 'input', 'native')
            packet['preprocessing_complete'] = False
            packet['optimization_admission'] = True
            record.write_text(json.dumps(packet))
            with self.assertRaisesRegex(ValueError, 'CANNOT_GRANT_ADMISSION'):
                module.verify_chain(paths, names, 'input', 'native')


if __name__ == '__main__':
    unittest.main()
