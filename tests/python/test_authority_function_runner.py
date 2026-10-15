"""The frozen runner must dispatch and persist authoritative, not float, state."""
import importlib.util
import json
import os
import stat
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

RUNNER = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"


def runner():
    spec = importlib.util.spec_from_file_location("authority_function_runner",RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AuthorityFunctionRunnerTests(unittest.TestCase):
    def test_terminal_attachment_is_durable_before_return(self):
        # Break caught: returning after buffered JSON output without syncing
        # the complete payload and its newly created directory entry.
        api = runner()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'native-run.json'
            synced = []
            real_sync = os.fsync
            def observe_sync(fd):
                mode = os.fstat(fd).st_mode
                synced.append('directory' if stat.S_ISDIR(mode) else 'file')
                if stat.S_ISREG(mode):
                    self.assertEqual(json.loads(path.read_text()),
                                     {'coefficient': '1/3', 'complete': True})
                real_sync(fd)
            with patch.object(api.os, 'fsync', side_effect=observe_sync):
                api._write(path, {'coefficient': '1/3', 'complete': True})
            self.assertEqual(synced, ['file', 'directory'])

    def test_terminal_attachment_sync_failure_is_not_success(self):
        # A material export must not continue if persistence failed. Keep the
        # file as diagnostic evidence, never overwrite it on an attempted retry.
        api = runner()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'native-run.json'
            with patch.object(api.os, 'fsync', side_effect=OSError('sync failed')):
                with self.assertRaisesRegex(OSError, 'sync failed'):
                    api._write(path, {'complete': True})
            self.assertTrue(path.exists())
            with self.assertRaises(FileExistsError):
                api._write(path, {'complete': False})
            self.assertEqual(json.loads(path.read_text()), {'complete': True})

    def test_terminal_attachment_directory_sync_failure_propagates(self):
        api = runner()
        real_sync = os.fsync
        def fail_directory(fd):
            if stat.S_ISDIR(os.fstat(fd).st_mode):
                raise OSError('directory sync failed')
            real_sync(fd)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'native-run.json'
            with patch.object(api.os, 'fsync', side_effect=fail_directory):
                with self.assertRaisesRegex(OSError, 'directory sync failed'):
                    api._write(path, {'complete': True})
            self.assertEqual(json.loads(path.read_text()), {'complete': True})

    def test_saved_complete_reinspection_rejects_diagnostic_and_wrong_frozen_request(self):
        # Break caught: trusting saved passed=true instead of replaying a frozen
        # full run. These intentionally incomplete files are NOT solver results.
        api=runner()
        self.assertTrue(callable(getattr(api,'reinspect_saved_authority_complete_material',None)),
                        'final audit must independently re-read complete saved material')
        request,objective,case=api.load_frozen_rotation('cuboid-offaxis-pairs-yoyo','YOYO',0)
        native=api.ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        evaluator=api.frozen_rotation_evaluator(exact_source_boundary=True,exact_material_boundary=True,
            rational_source_coordinates=True,compact_rational_source_coordinates=True)
        frozen=api.frozen_run_inputs(request,objective,case,0,native,evaluator)
        frozen.update(coefficient_policy='COEFFICIENT_MPFR1024_V1',coefficient_precision_bits=1024)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            def put(name,value): (root/name).write_text(json.dumps(value)+'\n')
            put('frozen-input.json',frozen)
            put('native-run.json',dict(schema='spin2014-authority-actual-point-v1',passed=True))
            with self.assertRaises(ValueError):api.reinspect_saved_authority_complete_material(root)
            frozen['fixed_material_density']=2
            put('frozen-input.json',frozen)
            put('native-run.json',dict(schema='spin2014-authority-restarts-v1',
                                      fixed_material_density=2,passed=True))
            with self.assertRaisesRegex(ValueError,'saved request differs'):
                api.reinspect_saved_authority_complete_material(root)

    def test_repeat_process_receipts_reject_cloned_units_and_truncated_jobs(self):
        # Receipt parsing only: these synthetic records are NOT real experiments.
        # The break caught is equality being upgraded without distinct complete
        # whole-workflow resource receipts, or with a changed product identity.
        import copy
        from make_it_stand.verification.source_checkpoint_store import _fingerprints
        api=runner()
        self.assertTrue(callable(getattr(api,'audit_authority_repeat_processes',None)),
                        'repeat admission needs independent bounded process receipts')
        native=api.ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        def put(path,value):path.write_text(json.dumps(value,sort_keys=True)+'\n')
        def write(root,unit):
            root.mkdir();(root/'output').mkdir()
            inputs={'runner':str(RUNNER),'native':str(native)}
            config=dict(schema='source-job-v1',unit=unit,boot_id='test-boot',
                identity_files=inputs,identity=_fingerprints(inputs),
                started_boottime_ns=1000,deadline_boottime_ns=3600000001000,
                limits=dict(maximum_wall_ns=3600000000000,maximum_memory_bytes=24*1024**3),
                argv=[sys.executable,'-B',str(RUNNER),'--case','cuboid-offaxis-pairs-yoyo',
                      '--objective','YOYO','--seed','0','--coefficient-policy',
                      'COEFFICIENT_MPFR1024_V1','--output',str(root/'output')])
            group='/user.slice/'+unit
            result=dict(schema='source-job-result-v1',unit=unit,status='PROCESS_EXITED',
                returncode=0,elapsed_ns=1000000000,service=dict(Result='success',
                SubState='exited',ExecMainCode='1',ExecMainStatus='0',ControlGroup=group,
                MemoryMax=str(24*1024**3),MemorySwapMax='0',MemorySwapPeak='0'))
            kernel=dict(boot_id='test-boot',cgroup=group,memory_max=24*1024**3,swap_max=0,
                guardian_pid=100,observed_boottime_ns=2000,
                absolute_deadline_boottime_ns=config['deadline_boottime_ns'])
            put(root/'config.json',config);put(root/'result.json',result)
            put(root/'complete.json',dict(result_sha256=api._sha(root/'result.json')))
            put(root/'kernel.json',kernel)
            put(root/'output/frozen-input.json',dict(case='cuboid-offaxis-pairs-yoyo',
                objective=dict(objective='YOYO'),seed=0,
                coefficient_policy='COEFFICIENT_MPFR1024_V1',coefficient_precision_bits=1024,
                native_path=str(native),native_sha256=api._sha(native)))
            return config,result,kernel
        with tempfile.TemporaryDirectory() as directory:
            first=Path(directory)/'first';repeat=Path(directory)/'repeat'
            write(first,'first.service');config,result,kernel=write(repeat,'repeat.service')
            checked=api.audit_authority_repeat_processes(first,repeat)
            self.assertTrue(checked['independent_bounded_jobs'])
            self.assertFalse(checked['sealed_artifact'])
            self.assertEqual(checked['scientific_material_admission'],'NOT_PERFORMED')
            bad=copy.deepcopy(result);bad['status']='RESOURCE_LIMIT_WALL'
            put(repeat/'result.json',bad)
            put(repeat/'complete.json',dict(result_sha256=api._sha(repeat/'result.json')))
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)
            put(repeat/'result.json',result)
            put(repeat/'complete.json',dict(result_sha256=api._sha(repeat/'result.json')))
            bad_config=copy.deepcopy(config);bad_config['unit']='first.service'
            bad_result=copy.deepcopy(result);bad_result['unit']='first.service'
            bad_result['service']['ControlGroup']='/user.slice/first.service'
            bad_kernel=copy.deepcopy(kernel);bad_kernel['cgroup']='/user.slice/first.service'
            put(repeat/'config.json',bad_config);put(repeat/'result.json',bad_result)
            put(repeat/'complete.json',dict(result_sha256=api._sha(repeat/'result.json')))
            put(repeat/'kernel.json',bad_kernel)
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)
            put(repeat/'config.json',config);put(repeat/'result.json',result)
            put(repeat/'kernel.json',kernel)
            put(repeat/'complete.json',dict(result_sha256=api._sha(repeat/'result.json')))
            # A reboot is not itself evidence of a different scientific
            # environment. Keep boot identity for process disambiguation,
            # not as an extra scientific completion requirement.
            other_config=copy.deepcopy(config);other_config['boot_id']='other-boot'
            other_kernel=copy.deepcopy(kernel);other_kernel['boot_id']='other-boot'
            put(repeat/'config.json',other_config);put(repeat/'kernel.json',other_kernel)
            self.assertTrue(api.audit_authority_repeat_processes(first,repeat)['independent_bounded_jobs'])
            put(repeat/'config.json',config);put(repeat/'kernel.json',kernel)
            frozen=json.loads((repeat/'output/frozen-input.json').read_text())
            bad=copy.deepcopy(frozen);bad['native_sha256']='0'*64
            put(repeat/'output/frozen-input.json',bad)
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)
            put(repeat/'output/frozen-input.json',frozen)
            bad=copy.deepcopy(config);bad['argv'][bad['argv'].index('--seed')+1]='1'
            put(repeat/'config.json',bad)
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)
            put(repeat/'config.json',config)
            bad=copy.deepcopy(config);bad['argv'][0]='/usr/bin/true'
            put(repeat/'config.json',bad)
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)
            put(repeat/'config.json',config)
            bad=copy.deepcopy(kernel);bad['swap_max']=False
            put(repeat/'kernel.json',bad)
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)
            put(repeat/'kernel.json',kernel)
            bad=copy.deepcopy(result);bad['service']['MemorySwapPeak']='1'
            put(repeat/'result.json',bad)
            put(repeat/'complete.json',dict(result_sha256=api._sha(repeat/'result.json')))
            with self.assertRaises(ValueError):api.audit_authority_repeat_processes(first,repeat)

    def test_event_attempt_agreement_rejects_missing_updates_and_wrong_trial(self):
        # Structural observation fixture only; never a material/solver pass.
        import copy
        api=runner()
        self.assertTrue(callable(getattr(api,'audit_authority_event_attempt_agreement',None)),
                        'returned attempts must agree with live scientific events')
        trial=dict(arc_length='1/2',reason='residual_decreased',diagnostic='',
                   raw_equation_residual='1/100',trial_identity='t'*64,
                   geometry_authority=['2/1'],dual_authority=['3/1'],accepted=True)
        attempt=dict(function_evaluations=2,accepted_updates=1,
            initial_geometry_authority=['1/1'],initial_dual_authority=['1/1'],
            initial_state_identity='i'*64,initial_material_identity='m'*64,
            initial_raw_residual='1/1',accepted_state_identities=['f'*64],trials=[trial],
            final_state_identity='f'*64,final_raw_residual='1/100',
            status='stationary_candidate',diagnostic='')
        def event(kind,index,**values):
            return dict(type=kind,restart=0,evaluation_index=index,**values)
        events=[event('attempt_start',0,geometry_authority=['1/1'],dual_authority=['1/1']),
            event('initial_evaluation_start',1),
            event('initial_evaluation_terminal',1,state_identity='i'*64,
                  material_identity='m'*64,raw_residual='1/1',status='',diagnostic=''),
            event('trial_evaluation_start',2,geometry_authority=['2/1'],
                  dual_authority=['3/1'],arc_length='1/2',state_identity='i'*64),
            event('trial_evaluation_terminal',2,geometry_authority=['2/1'],
                  dual_authority=['3/1'],arc_length='1/2',status='residual_decreased',
                  diagnostic='',raw_residual='1/100',trial_identity='t'*64,accepted=True,
                  state_identity='i'*64),
            event('state_accepted',2,state_identity='f'*64,raw_residual='1/100',
                  geometry_authority=['2/1'],dual_authority=['3/1'],accepted=True,
                  material_identity='z'*64),
            event('attempt_terminal',2,status='stationary_candidate',diagnostic='',
                  state_identity='f'*64,raw_residual='1/100'),
            dict(type='run_terminal',restart=1,evaluation_index=2,status='stationary_candidate')]
        for index,item in enumerate(events): item['sequence']=index
        payload=dict(attempts=[attempt],events=events,function_evaluations=2,status='stationary_candidate',
                     stationary_point=dict(material_scientific_hash='z'*64))
        checked=api.audit_authority_event_attempt_agreement(payload)
        self.assertTrue(checked['events_match_attempts'])
        self.assertFalse(checked['sealed_artifact'])
        bad=copy.deepcopy(payload); bad['events'].pop(5)
        for index,item in enumerate(bad['events']): item['sequence']=index
        with self.assertRaises(ValueError): api.audit_authority_event_attempt_agreement(bad)
        bad=copy.deepcopy(payload);bad['events'][4]['raw_residual']='1/10'
        with self.assertRaises(ValueError): api.audit_authority_event_attempt_agreement(bad)
        bad=copy.deepcopy(payload);bad['events'][-1]['evaluation_index']=1
        with self.assertRaises(ValueError): api.audit_authority_event_attempt_agreement(bad)
        for index,key,value in ((3,'state_identity','other-base'),
                (4,'state_identity','other-base'),
                (2,'status','initial_geometry_or_derivative_failure'),
                (5,'material_identity','other-material')):
            bad=copy.deepcopy(payload);bad['events'][index][key]=value
            with self.subTest(index=index,key=key),self.assertRaises(ValueError):
                api.audit_authority_event_attempt_agreement(bad)
        # A real construction failure consumes one evaluation, then a fresh
        # attempt may succeed; do not reinterpret failed starts as accepted work.
        retried=copy.deepcopy(payload)
        failed=dict(function_evaluations=1,accepted_updates=0,trials=[],
            accepted_state_identities=[],initial_geometry_authority=['0/1'],
            initial_dual_authority=['0/1'],initial_state_identity='',
            initial_material_identity='',initial_raw_residual='',
            final_state_identity='',final_raw_residual='',
            status='initial_geometry_or_derivative_failure',diagnostic='construction failed')
        failed_events=[event('attempt_start',0,geometry_authority=['0/1'],dual_authority=['0/1']),
            event('initial_evaluation_start',1),
            event('initial_evaluation_terminal',1,state_identity='',material_identity='',
                  raw_residual='',status=failed['status'],diagnostic=failed['diagnostic']),
            event('attempt_terminal',1,state_identity='',raw_residual='',
                  status=failed['status'],diagnostic=failed['diagnostic'])]
        for item in retried['events']:item['restart']+=1
        retried['attempts'].insert(0,failed)
        retried['events']=failed_events+retried['events']
        retried['function_evaluations']=3
        retried['events'][-1]['evaluation_index']=3
        for index,item in enumerate(retried['events']):item['sequence']=index
        self.assertTrue(api.audit_authority_event_attempt_agreement(retried)['events_match_attempts'])

    def test_material_witness_failure_preserves_diagnostic_without_exclusive_rewrite(self):
        # Persistence boundary only, not a scientific pass: substitute the seed
        # reconstruction, but execute the real malformed-moment check and file IO.
        api=runner()
        self.assertTrue(callable(getattr(api,'record_authority_seed_audit',None)),
                        'runner needs a failure-safe seed audit recording boundary')
        payload=dict(attempts=[dict(accepted_updates=1,initial_total_moments=[])],
                     stationary_attempt_index=0,stationary_point=dict(total_moments_target=[]))
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory); summary={}
            with patch.object(api,'audit_authority_seed_ledger',return_value=dict(
                    full_run_candidate=True,sealed_artifact=False)):
                api.record_authority_seed_audit(None,None,payload,output,summary)
            diagnostic=json.loads((output/'independent-seed-ledger-audit.json').read_text())
            self.assertFalse(summary['seed_and_restart_ledger_verified'])
            self.assertEqual(diagnostic['exception_type'],'ValueError')
            self.assertIn('moment',diagnostic['diagnostic'])
            self.assertFalse(diagnostic['sealed_artifact'])
            self.assertFalse((output/'material-change-moment-witness.json').exists())

    def test_repeat_comparison_binds_saved_authority_and_preserves_failed_attempts(self):
        # Break caught: comparing only labels, ignoring exact coefficients,
        # losing failed attempts, or comparing an edited scientific projection.
        api=runner()
        def put(path,value):
            path.write_text(json.dumps(value,sort_keys=True,allow_nan=False)+'\n')
        def write(root,*,seconds):
            identity=dict(request_identity='a'*64,context_identity='b'*64,objective_identity='c'*64)
            events=[dict(identity,type='attempt_start',sequence=0,
                         geometry_authority=['1/100000000000000000000000000000001']),
                    dict(identity,type='run_terminal',sequence=1,status='restart_budget_exhausted')]
            payload=dict(schema='spin2014-authority-restarts-v1',
                request_identity='a'*64,context_identity='b'*64,
                protocol_scope='CONFIGURED_FULL_ATTEMPT_BUDGET',events=events,
                status='restart_budget_exhausted',attempts=[dict(status='initial_construction_failed')],
                evaluator_wall_seconds=seconds)
            put(root/'native-run.json',payload)
            put(root/'scientific-output.json',api._scientific(payload))
            put(root/'frozen-input.json',dict(coefficient_policy='COEFFICIENT_MPFR1024_V1',
                coefficient_precision_bits=1024,seed=0,fixed_material_density=2))
            with (root/'solver-events.jsonl').open('w') as stream:
                for i,event in enumerate(events):
                    stream.write(json.dumps(dict(schema='spin2014-solver-event-journal-v1',sequence=i,event=event))+'\n')
            return payload
        with tempfile.TemporaryDirectory() as directory:
            first=Path(directory)/'first';repeat=Path(directory)/'repeat'
            first.mkdir();repeat.mkdir()
            write(first,seconds=1);payload=write(repeat,seconds=2)
            checked=api.compare_authority_scientific_outputs(first,repeat)
            self.assertTrue(checked['scientific_outputs_identical'])
            self.assertFalse(checked['sealed_artifact'])
            self.assertEqual(checked['final_numeric_material_admission'],'NOT_PERFORMED')
            with self.assertRaises(ValueError):
                api.compare_authority_scientific_outputs(first,first)
            # A smaller denominator difference must not disappear through float views.
            payload['events'][0]['geometry_authority']=['1/100000000000000000000000000000002']
            put(repeat/'native-run.json',payload)
            with self.assertRaises(ValueError):
                api.compare_authority_scientific_outputs(first,repeat)
            put(repeat/'scientific-output.json',api._scientific(payload))
            with (repeat/'solver-events.jsonl').open('w') as stream:
                for i,event in enumerate(payload['events']):
                    stream.write(json.dumps(dict(schema='spin2014-solver-event-journal-v1',sequence=i,event=event))+'\n')
            checked=api.compare_authority_scientific_outputs(first,repeat)
            self.assertFalse(checked['scientific_outputs_identical'])
            self.assertNotEqual(checked['scientific_sha256'][0],checked['scientific_sha256'][1])
            payload=write(repeat,seconds=2)
            payload['attempts']=[]
            put(repeat/'native-run.json',payload)
            put(repeat/'scientific-output.json',api._scientific(payload))
            checked=api.compare_authority_scientific_outputs(first,repeat)
            self.assertFalse(checked['scientific_outputs_identical'])
            write(repeat,seconds=2)
            frozen=json.loads((repeat/'frozen-input.json').read_text())
            frozen['fixed_material_density']=3
            put(repeat/'frozen-input.json',frozen)
            checked=api.compare_authority_scientific_outputs(first,repeat)
            self.assertFalse(checked['frozen_inputs_identical'])

    def test_complete_export_rejects_returned_density_not_frozen_request(self):
        # Break caught: rebuilding an initial material with a returned, rather
        # than frozen, fixed-shell density authenticates the wrong problem.
        api=runner()
        request,objective,_=api.load_frozen_rotation('cuboid-offaxis-pairs-yoyo','YOYO',0)
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            payload=dict(status='stationary_candidate',
                protocol_scope='CONFIGURED_FULL_ATTEMPT_BUDGET',fixed_material_density=2)
            (output/'native-run.json').write_text(json.dumps(payload))
            (output/'frozen-input.json').write_text(json.dumps(dict(fixed_material_density=1)))
            with self.assertRaises(Exception) as caught:
                api.export_authority_complete_material(request,objective,None,payload,output)
            self.assertIsInstance(caught.exception,ValueError)
            self.assertIn('fixed density differs from frozen request',str(caught.exception))
            self.assertFalse((output/'authority-initial-point.json').exists())

    def test_authority_audit_files_are_in_run_evidence_outputs(self):
        api=runner()
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            names=('independent-event-journal-audit.json','independent-seed-ledger-audit.json',
                   'material-change-moment-witness.json')
            for name in names:
                (output/name).write_text('{}')
            self.assertEqual(set(api.authority_audit_outputs(output)),{output/name for name in names})

    def test_complete_material_export_rejects_failed_run_before_output(self):
        api=runner()
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'not-created'
            with self.assertRaisesRegex(ValueError,'stationary full run required'):
                api.export_authority_complete_material(None,None,None,
                    dict(status='restart_budget_exhausted'),output)
            self.assertFalse(output.exists())
    def test_seed_audit_reconstructs_starts_and_rejects_changed_start(self):
        # Structural fixture only: no material or optimization success claimed.
        import copy
        from make_it_stand.spin import _native_request
        api=runner()
        request,objective,case=api.load_frozen_rotation('cuboid-offaxis-pairs-yoyo','YOYO',0)
        native,nr=_native_request(request)
        starts=native.spin2014_authority_initial_sequence(nr.mode,nr.settings)
        attempts=[dict(initial_geometry_authority=s['geometry_authority'],
            initial_dual_authority=s['dual_authority'],function_evaluations=1,
            initial_state_identity='',trials=[],accepted_updates=0,
            accepted_state_identities=[],step_identities=[],status='initial_construction_failed') for s in starts]
        payload=dict(schema='spin2014-authority-restarts-v1',request_identity='a'*64,
            context_identity='b'*64,protocol_scope='CONFIGURED_FULL_ATTEMPT_BUDGET',
            attempts=attempts,stationary_attempt_index=None,stationary_point=None,
            status='restart_budget_exhausted',function_evaluations=len(attempts))
        checked=api.audit_authority_seed_ledger(request,objective,payload)
        self.assertFalse(checked['full_run_candidate'])
        self.assertFalse(checked['sealed_artifact'])
        bad=copy.deepcopy(payload)
        original=bad['attempts'][0]['initial_geometry_authority'][0]
        bad['attempts'][0]['initial_geometry_authority'][0]='1/1' if original=='0/1' else '0/1'
        with self.assertRaisesRegex(ValueError,'initial sequence mismatch'):
            api.audit_authority_seed_ledger(request,objective,bad)

    def test_policy_cannot_relabel_saved_or_static_run(self):
        # Break caught: using the MPFR policy label on a historical float export.
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/"not-created"
            for args in (("--saved-static-run","missing.json"),
                         ("--saved-rotational-run","missing.json"),
                         ("--case","static-faceted-shell-lowering","--objective","STATIC_LOWERING","--seed","0")):
                result = subprocess.run([sys.executable,"-B",str(RUNNER),*args,
                    "--coefficient-policy","COEFFICIENT_MPFR1024_V1","--output",str(output)],
                    capture_output=True,text=True,timeout=20)
                self.assertEqual(result.returncode,2)
                self.assertIn("coefficient policy requires a fresh rotational run",result.stderr)
                self.assertFalse(output.exists())

    def test_payload_keeps_rational_authority_and_diagnostic_scope(self):
        # Break caught: passing authority through asdict/float material-run conversion.
        value = dict(schema="spin2014-authority-restarts-v1",
            scope="NUMERICAL_CANDIDATE_NOT_FINAL_SEAL",attempts=[],
            stationary_point={"geometry_authority":["1/1000000000000000000000000000000000000000000000000000001"]})
        payload = runner().solver_run_payload(value)
        self.assertEqual(json.loads(json.dumps(payload)),value)
        self.assertNotIn("material_export_complete",payload)

    def test_interrupted_run_does_not_report_zero_executed_evaluations(self):
        # Break caught: defaulting a missing native return to zero actual cost.
        api = runner()
        partial = api.authority_run_accounting(dict(schema="spin2014-quadric-admission-failure-v1"))
        self.assertIsNone(partial["function_evaluations"])
        self.assertIsNone(partial["accepted_steps"])
        self.assertFalse(partial["native_return_complete"])
        complete = api.authority_run_accounting(dict(schema="spin2014-authority-restarts-v1",
            function_evaluations=3,attempts=[dict(function_evaluations=1,accepted_updates=0),
            dict(function_evaluations=2,accepted_updates=1)],events=[dict(type="run_terminal")]))
        self.assertEqual(complete["function_evaluations"],3)
        self.assertEqual(complete["restart_function_evaluations"],[1,2])
        self.assertEqual(complete["accepted_steps"],1)
        self.assertTrue(complete["native_return_complete"])

    def test_failed_screens_and_budget_terminal_are_not_ok_events(self):
        api = runner()
        for status in ("initial_multiplier_nonpositive","evaluation_budget_exhausted",
                       "restart_budget_exhausted","linear_solve_or_direction_failed"):
            self.assertEqual(api.authority_event_outcome(dict(type="attempt_terminal",status=status)),"failed")
        self.assertEqual(api.authority_event_outcome(dict(type="trial_evaluation_terminal",
            status="no_residual_decrease")),"failed")
        self.assertEqual(api.authority_event_outcome(dict(type="trial_evaluation_terminal",
            status="topology_transition_relinearized")),"ok")


if __name__ == "__main__":
    unittest.main()
