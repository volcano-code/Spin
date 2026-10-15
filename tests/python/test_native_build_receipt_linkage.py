"""Individually valid observations from two builds must not be spliced."""
import json
import shutil
import unittest

import test_native_build_receipt as fixture


class NativeBuildReceiptLinkageTests(unittest.TestCase):
    def setUp(self):
        self.f = fixture.NativeBuildReceiptTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        for path in (self.f.config, self.f.result):
            data = json.loads(path.read_text())
            data['unit'] = 'build-a.service'
            fixture.save(path, data)
        fixture.save(self.f.complete, {'result_sha256': fixture.sha(self.f.result)})
        self.f.refresh()

    def test_same_job_and_external_reference_remain_valid(self):
        explicit = self.f.root / 'external-receipt.json'
        self.f.receipt.rename(explicit)
        self.assertIn(self.f.cache, self.f.resolve(explicit))

    def test_valid_second_build_summary_cannot_replace_first_build_product(self):
        f = self.f
        other_job = f.root / 'job-b'
        shutil.copytree(f.job, other_job)
        other_build = other_job / 'output/bundle'
        other_module = other_build / 'python/native.so'
        other_summary = other_job / 'output/summary.json'
        fixture.save(other_summary, {'status': 'OBSERVED_COMPILE_ARCHIVE_LINK_LOAD',
            'module': {'path': str(other_module), 'sha256': fixture.sha(other_module)}})
        receipt = json.loads(f.receipt.read_text())
        receipt['files']['summary'] = {'path': str(other_summary), 'sha256': fixture.sha(other_summary)}
        mixed = other_build / 'build-receipt.json'
        fixture.save(mixed, receipt)
        with self.assertRaisesRegex(ValueError, 'same build|layout|job'):
            fixture.run_evidence.native_build_configuration_files(other_build, other_module)

    def test_different_job_unit_is_rejected_with_consistent_file_hashes(self):
        result = json.loads(self.f.result.read_text())
        result['unit'] = 'build-b.service'
        fixture.save(self.f.result, result)
        fixture.save(self.f.complete, {'result_sha256': fixture.sha(self.f.result)})
        self.f.refresh()
        with self.assertRaisesRegex(ValueError, 'job|unit'):
            self.f.resolve()

    def test_missing_job_identity_is_rejected(self):
        config = json.loads(self.f.config.read_text())
        del config['unit']
        fixture.save(self.f.config, config)
        self.f.refresh()
        with self.assertRaisesRegex(ValueError, 'job|unit'):
            self.f.resolve()

