"""Receipt binding must not turn an unrelated cache/library into build evidence."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from make_it_stand.verification import run_evidence


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding='utf-8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class NativeBuildReceiptTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.build = self.root/'job/output/bundle'
        self.module = self.build/'python/native.so'
        self.module.parent.mkdir(parents=True)
        self.module.write_bytes(b'fixture binary, not executable')
        self.cache = self.root/'configured/CMakeCache.txt'
        self.cache.parent.mkdir()
        self.cache.write_text('CMAKE_BUILD_TYPE:STRING=Release\n')
        self.job = self.root/'job'
        self.config = self.job/'config.json'
        self.result = self.job/'result.json'
        self.complete = self.job/'complete.json'
        self.summary = self.job/'output/summary.json'
        self.receipt = self.build/'build-receipt.json'
        save(self.config, {'schema':'source-job-v1', 'unit':'fixture-build.service', 'identity':{'cache':{
            'path':str(self.cache), 'sha256':sha(self.cache), 'size':self.cache.stat().st_size}}})
        save(self.result, {'schema':'source-job-result-v1','unit':'fixture-build.service',
            'status':'PROCESS_EXITED','returncode':0})
        save(self.complete, {'result_sha256':sha(self.result)})
        save(self.summary, {'status':'OBSERVED_COMPILE_ARCHIVE_LINK_LOAD',
            'module':{'path':str(self.module),'sha256':sha(self.module)}})
        self.refresh()

    def refresh(self):
        save(self.receipt, {'schema':'native-build-receipt-reference-v1',
            'files':{key:{'path':str(path),'sha256':sha(path)} for key,path in (
                ('config',self.config),('result',self.result),('complete',self.complete),('summary',self.summary))}})

    def resolve(self, receipt=None):
        self.assertTrue(hasattr(run_evidence, 'native_build_configuration_files'),
                        'actual runtime cannot bind observed build receipt')
        return run_evidence.native_build_configuration_files(self.build, self.module, receipt=receipt)

    def test_real_cache_is_selected_and_every_receipt_input_is_returned(self):
        self.assertEqual(set(self.resolve()), {self.receipt,self.config,self.result,self.complete,self.summary,self.cache})
        self.assertFalse((self.build/'CMakeCache.txt').exists())

    def test_legacy_cache_remains_supported_without_receipt(self):
        self.receipt.unlink()
        local = self.build/'CMakeCache.txt'
        local.write_text('legacy configuration')
        self.assertEqual(self.resolve(), [local])

    def test_explicit_receipt_outside_bundle_is_supported(self):
        explicit=self.root/'explicit.json'
        self.receipt.rename(explicit)
        self.assertIn(explicit,self.resolve(explicit))

    def test_changed_cache_is_rejected(self):
        self.cache.write_text('different configuration')
        with self.assertRaisesRegex(ValueError,'cache|identity'):self.resolve()

    def test_changed_loaded_binary_is_rejected(self):
        self.module.write_bytes(b'wrong executable')
        with self.assertRaisesRegex(ValueError,'module|identity'):self.resolve()

    def test_unrelated_module_path_is_rejected_even_with_identical_bytes(self):
        other=self.root/'other.so';other.write_bytes(self.module.read_bytes())
        save(self.summary,{'status':'OBSERVED_COMPILE_ARCHIVE_LINK_LOAD','module':{'path':str(other),'sha256':sha(other)}})
        self.refresh()
        with self.assertRaisesRegex(ValueError,'module'):self.resolve()

    def test_failed_build_cannot_be_granted_by_a_refreshed_reference(self):
        save(self.result,{'schema':'source-job-result-v1','status':'PROCESS_FAILED','returncode':1})
        save(self.complete,{'result_sha256':sha(self.result)})
        self.refresh()
        with self.assertRaisesRegex(ValueError,'successful'):self.resolve()

    def test_incomplete_or_mismatched_completion_is_rejected(self):
        save(self.complete,{'result_sha256':'0'*64});self.refresh()
        with self.assertRaisesRegex(ValueError,'completion'):self.resolve()

    def test_mutated_receipt_dependency_and_symlink_are_rejected(self):
        self.config.write_text('{}')
        with self.assertRaisesRegex(ValueError,'identity'):self.resolve()
        self.refresh()
        original=self.config.with_suffix('.saved');self.config.rename(original)
        self.config.symlink_to(original)
        with self.assertRaisesRegex(ValueError,'symlink'):self.resolve()

    def test_ambiguous_cache_is_rejected(self):
        data=json.loads(self.config.read_text())
        data['identity']['second']=data['identity']['cache']
        save(self.config,data);self.refresh()
        with self.assertRaisesRegex(ValueError,'exactly one'):self.resolve()

    def test_invalid_receipt_does_not_fall_back_to_a_local_cache(self):
        (self.build/'CMakeCache.txt').write_text('unbound local cache')
        save(self.receipt,{'schema':'bad'})
        with self.assertRaisesRegex(ValueError,'receipt'):self.resolve()
