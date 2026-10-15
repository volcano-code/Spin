"""Split-job accounting tests; synthetic packets are NOT derivative evidence."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import subprocess
import sys
import os

from make_it_stand.verification import quadric_point_jobs as jobs


class QuadricPointJobs(unittest.TestCase):
    def test_freezer_refuses_unbound_original_registration_before_native_or_output(self):
        script=Path(__file__).resolve().parents[2]/'experiments/spin2014_objectives_quadric_v1/freeze_production_points.py'
        with tempfile.TemporaryDirectory() as temp:
            original=Path(temp)/'original.json';original.write_text('{}')
            output=Path(temp)/'must-not-exist.json'
            result=subprocess.run([sys.executable,'-B',str(script),'--registration',str(original),
                '--registration-sha256','a'*64,'--output',str(output)],
                capture_output=True,text=True,env={**os.environ,'PYTHONPATH':str(script.parents[2]/'python')})
            self.assertNotEqual(result.returncode,0)
            self.assertIn('original registration identity mismatch',result.stderr)
            self.assertFalse(output.exists())

    def test_worker_without_valid_resource_envelope_refuses_before_input_or_output(self):
        # Admission guard must run BEFORE missing input read / Boolean / writes.
        script=Path(__file__).resolve().parents[2]/'experiments/spin2014_objectives_quadric_v1/production_point_job.py'
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/'must-not-be-created'
            result=subprocess.run([sys.executable,'-B',str(script),'--manifest',str(Path(temp)/'missing'),
                '--manifest-identity','a'*64,'--index','0','--output',str(output),
                '--job-envelope',str(Path(temp)/'kernel.json')],
                capture_output=True,text=True,env={**os.environ,'PYTHONPATH':str(script.parents[2]/'python')})
            self.assertNotEqual(result.returncode,0)
            self.assertRegex(result.stderr,'24 GiB / swap 0|SOURCE_JOB_RECEIPT_UNAVAILABLE')
            self.assertFalse(output.exists())

    def test_aggregator_without_valid_resource_envelope_cannot_create_pass_artifact(self):
        script=Path(__file__).resolve().parents[2]/'experiments/spin2014_objectives_quadric_v1/aggregate_production_points.py'
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/'must-not-be-created'
            result=subprocess.run([sys.executable,'-B',str(script),'--manifest',str(Path(temp)/'missing'),
                '--manifest-identity','a'*64,'--jobs',temp,'--output',str(output),
                '--job-envelope',str(Path(temp)/'kernel.json')],
                capture_output=True,text=True,env={**os.environ,'PYTHONPATH':str(script.parents[2]/'python')})
            self.assertNotEqual(result.returncode,0)
            self.assertRegex(result.stderr,'24 GiB / swap 0|SOURCE_JOB_RECEIPT_UNAVAILABLE')
            self.assertFalse(output.exists())

    def manifest(self):
        return jobs.freeze_point_manifest([.25]*10,[.001]*10,
            sources={name:'a'*64 for name in ('input','native','runner','checker')},
            mapping_settings={'coefficient_policy':'COEFFICIENT_MPFR1024_V1',
                              'normal_level':2,'precision_bits':1024})

    def test_manifest_retains_all_shared_nodes_without_repeating_base(self):
        manifest=self.manifest()
        self.assertEqual(len(manifest['points']),81)
        self.assertEqual(len(manifest['rows']),30)
        self.assertEqual(manifest['rows'][0]['indices'],[1,3,0,4,2])
        self.assertEqual(manifest['rows'][2]['indices'],[5,7,0,8,6])
        self.assertEqual(manifest['points'][1]['raw_coefficients'][0],.249)
        self.assertEqual(manifest['points'][8]['raw_coefficients'][0],.250125)
        self.assertFalse(manifest['scientific_gate_granted'])

    def test_mutated_manifest_and_missing_binary_identity_rejected(self):
        manifest=self.manifest()
        manifest['points'][1]['raw_coefficients'][0]=.5
        with self.assertRaises(ValueError): jobs.validate_point_manifest(manifest)
        with self.assertRaises(ValueError):
            jobs.freeze_point_manifest([.25]*10,[.001]*10,sources={'input':'a'*64},
                                       mapping_settings={})

    def test_source_change_rejected_before_point_execution(self):
        with tempfile.TemporaryDirectory() as temp:
            paths={name:Path(temp)/name for name in ('input','native','runner','checker')}
            for name,path in paths.items(): path.write_text(name)
            manifest=jobs.freeze_point_manifest([.25]*10,[.001]*10,
                sources={name:jobs.file_sha256(path) for name,path in paths.items()},
                mapping_settings={'coefficient_policy':'COEFFICIENT_MPFR1024_V1',
                                  'normal_level':2,'precision_bits':1024})
            jobs.verify_point_sources(manifest,paths)
            paths['native'].write_text('different loaded build')
            with self.assertRaises(ValueError): jobs.verify_point_sources(manifest,paths)
            with self.assertRaises(ValueError): jobs.verify_point_sources(manifest,{'input':paths['input']})

    def write_job(self,root,manifest,index, *, status='COMPLETE'):
        directory=root/str(index);directory.mkdir()
        packet={'index':index,'raw_coefficients':manifest['points'][index]['raw_coefficients'],
                'evaluation':{'scope':'SYNTHETIC_ACCOUNTING_TEST_NOT_PRODUCTION'}}
        path=directory/'point.json'
        path.write_text(json.dumps(packet))
        independent=directory/'independent.json'
        independent.write_text(json.dumps({'scope':'SYNTHETIC_ACCOUNTING_ONLY'}))
        receipt={'schema':'quadric-point-job-terminal-v1','manifest_sha256':manifest['identity'],
                 'index':index,'status':status,'actual_material_evaluations':1,
                 'packet_sha256':jobs.file_sha256(path),
                 'independent_sha256':jobs.file_sha256(independent)}
        (directory/'terminal.json').write_text(json.dumps(receipt))
        return directory

    def test_saved_packets_loaded_in_manifest_order_not_job_finish_order(self):
        manifest=self.manifest()
        with tempfile.TemporaryDirectory() as temp:
            directories=[self.write_job(Path(temp),manifest,i) for i in range(81)]
            packets=jobs.load_complete_point_jobs(manifest,list(reversed(directories)))
            self.assertEqual([p['index'] for p in packets],list(range(81)))
            self.assertEqual(packets[80]['raw_coefficients'][-1],.250125)
            # This loader grants only complete accounting, not geometry or FD success.
            self.assertNotIn('passed',packets[0])

    def test_stream_reopens_saved_material_and_rejects_modified_independent_attachment(self):
        # Catch eager loading all 81 large packets, or ignoring checker attachments.
        manifest=self.manifest()
        with tempfile.TemporaryDirectory() as temp:
            directories=[self.write_job(Path(temp),manifest,i) for i in range(81)]
            stream=jobs.iter_complete_point_jobs(manifest,directories)
            self.assertEqual(next(stream)['index'],0)
            (directories[1]/'independent.json').write_text('{}')
            with self.assertRaises(ValueError): next(stream)

    def test_missing_duplicate_resource_cut_and_stale_build_rejected(self):
        manifest=self.manifest()
        with tempfile.TemporaryDirectory() as temp:
            directories=[self.write_job(Path(temp),manifest,i) for i in range(81)]
            for bad in (directories[:-1],directories[:-1]+[directories[0]]):
                with self.assertRaises(ValueError): jobs.load_complete_point_jobs(manifest,bad)
            terminal=directories[-1]/'terminal.json'
            original=json.loads(terminal.read_text())
            for key,value in (('status','RESOURCE_LIMIT'),('manifest_sha256','b'*64),
                              ('actual_material_evaluations',0)):
                changed=copy.deepcopy(original);changed[key]=value
                terminal.write_text(json.dumps(changed))
                with self.assertRaises(ValueError): jobs.load_complete_point_jobs(manifest,directories)
            terminal.write_text(json.dumps(original))
            (directories[0]/'point.json').write_text('{}')
            with self.assertRaises(ValueError): jobs.load_complete_point_jobs(manifest,directories)


if __name__=='__main__': unittest.main()
