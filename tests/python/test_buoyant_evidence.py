"""Strict buoyant evidence admission must use original bounded native runs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from test_faster_cherchi_table_dispatch import strict_audit
from test_buoyant_table_dispatch import buoyant_case


class BuoyantEvidenceTests(unittest.TestCase):
    def test_saved_ledger_binds_request_labels_and_all_stage_events(self):
        from make_it_stand.verification.buoyant_application import buoyant_case_plan
        from make_it_stand.verification.buoyant_evidence import _validate_saved_ledger
        inputs=[buoyant_case()]
        case=buoyant_case_plan(inputs[0],1,1,[0])
        case['asset_class']=inputs[0].get('asset_class')
        # Ledger-only fixture, never a scientific success or material witness.
        stage={'stage':'A','status':'FAILED','trace':[{'evaluation':0,'failure':'test'}],
               'artifact_scientific_digest':None}
        case['runs'][0]['stages']=[stage]
        events=[{'iteration':0,'status':'failed','data':{'stage':'A','trace':stage['trace'],
                                                      'material_digest':None}}]
        self.assertEqual(_validate_saved_ledger(inputs,[case],events),events)
        for change in ('input','seed','missing_run','missing_failure','material','asset'):
            altered=copy.deepcopy(case)
            if change=='input':altered['buoyant_problem']['design']['density']=2.
            elif change=='seed':altered['runs'][0]['seed']=1
            elif change=='missing_run':altered['runs']=[]
            elif change=='missing_failure':altered['runs'][0]['stages'][0]['trace']=[]
            elif change=='material':altered['runs'][0]['stages'][0]['artifact_scientific_digest']='wrong'
            else:altered['asset_class']='AUTHOR_ORIGINAL'
            with self.subTest(change=change),self.assertRaises(ValueError):
                _validate_saved_ledger(inputs,[altered],events)
        with self.assertRaises(ValueError):_validate_saved_ledger(inputs,[case],[])

    def test_source_job_receipt_checks_real_files_and_resource_links(self):
        from make_it_stand.verification.buoyant_evidence import _job_receipt
        from make_it_stand.verification.run_evidence import _record
        from make_it_stand.verification.source_checkpoint_store import _fingerprints
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);out=root/'outputs';out.mkdir()
            execution=out/'evidence.json';execution.write_text('receipt fixture, not a V2 result')
            app=out/'application.json';app.write_text('{}')
            files={name:root/name for name in ('worker.py','input.json','native.so')}
            for path in files.values():path.write_text('frozen fixture')
            identity=_fingerprints(files)
            request={'asset_files':[str(files['input.json'])],'runner_files':[str(files['worker.py'])]}
            obs={'files':{'assets':[identity['input.json']], 'runners':[identity['worker.py']]},
                 'runtime':{'native_modules':{'_make_it_stand':identity['native.so']}}}
            config={'schema':'source-job-v1','unit':'test.service','boot_id':'boot',
                    'limits':{'maximum_wall_ns':100,'maximum_memory_bytes':1024,
                              'maximum_disk_bytes':1000,'minimum_start_free_bytes':100*1024**3,
                              'minimum_free_bytes':20*1024**3},
                    'identity_files':{k:str(p) for k,p in files.items()},'identity':identity,
                    'argv':[str(files['worker.py']),'worker',str(files['input.json']),str(out)]}
            kernel={'boot_id':'boot','memory_max':1024,'swap_max':0,'cgroup':'/test/test.service'}
            result={'schema':'source-job-result-v1','status':'PROCESS_EXITED','returncode':0,
                    'unit':'test.service','elapsed_ns':90,'accounted_disk_bytes':100,
                    'service':{'Result':'success','ExecMainStatus':'0','MemoryMax':'1024',
                               'MemorySwapMax':'0','MemorySwapPeak':'0','MemoryPeak':'512',
                               'ControlGroup':'/test/test.service'}}
            def save(c,k,r):
                for name,value in (('config',c),('kernel',k),('result',r)):
                    (root/(name+'.json')).write_text(json.dumps(value))
                (root/'complete.json').write_text(json.dumps({'result_sha256':_record(root/'result.json')['sha256']}))
            save(config,kernel,result)
            self.assertEqual(_job_receipt(root,request,obs,execution,app)['unit'],'test.service')
            for change in ('failed','timeout','swap','peak','memory_limit','boot','output','unpinned','disk'):
                c,k,r=copy.deepcopy((config,kernel,result))
                if change=='failed':r['returncode']=1
                elif change=='timeout':r['elapsed_ns']=101
                elif change=='swap':r['service']['MemorySwapPeak']='1'
                elif change=='peak':r['service']['MemoryPeak']='1025'
                elif change=='memory_limit':c['limits']['maximum_memory_bytes']=2048
                elif change=='boot':k['boot_id']='other'
                elif change=='output':c['argv'][-1]=str(root/'copied-output')
                elif change=='unpinned':c['identity'].pop('native.so');c['identity_files'].pop('native.so')
                else:r['accounted_disk_bytes']=1001
                save(c,k,r)
                with self.subTest(change=change),self.assertRaises(ValueError):
                    _job_receipt(root,request,obs,execution,app)
            save(config,kernel,result)
            files['native.so'].write_text('changed binary')
            with self.assertRaisesRegex(ValueError,'STALE_RESULT'):
                _job_receipt(root,request,obs,execution,app)
            files['native.so'].write_text('frozen fixture')
            (root/'complete.json').write_text('{"result_sha256":"wrong"}')
            with self.assertRaisesRegex(ValueError,'completion digest'):
                _job_receipt(root,request,obs,execution,app)

    def test_registry_requires_two_independent_checks(self):
        with self.assertRaisesRegex(ValueError, 'two independent'):
            strict_audit()._supporting_checks({'faster.buoyant_saved_material': []})

    def test_existing_material_record_without_execution_is_not_admitted(self):
        from make_it_stand.verification import buoyant_evidence
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            app=root/'application.json';app.write_text('{"cases":[]}')
            with self.assertRaisesRegex(ValueError,'original execution'):
                buoyant_evidence.audit_buoyant_execution(application=app,bundle_root=root,
                    execution_evidence=root/'missing.json',source_job_root=root)

    def test_repeat_identity_rejects_copy_even_with_distinct_paths(self):
        from make_it_stand.verification import buoyant_evidence
        record={'execution_sha256':'a'*64,'boot_id':'boot','unit':'unit1'}
        with self.assertRaisesRegex(ValueError,'same native run'):
            buoyant_evidence.assert_independent_buoyant_runs([record,dict(record)])
        other={**record,'execution_sha256':'b'*64}
        with self.assertRaisesRegex(ValueError,'same native run'):
            buoyant_evidence.assert_independent_buoyant_runs([record,other])


if __name__=='__main__':unittest.main()
