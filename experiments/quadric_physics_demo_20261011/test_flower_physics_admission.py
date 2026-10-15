import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from flower_run_physics import admit


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


class FlowerPhysicsAdmissionTests(unittest.TestCase):
    def test_real_flexible_result_merges_matching_scene_timestep_once(self):
        import flower_run_physics as api
        from flexible_bearing import support_protocol,simulate
        # A two-step physical kernel regression, not the full scene matrix.
        result=simulate(.005,[0.,0.,0.],np.diag([1e-6,1e-6,2e-6]),
                        protocol=support_protocol(.005,.002),dt=.0005,duration=.001,hz=35.,phase=0.)
        metadata={k:v for k,v in result.items() if k!='records'}
        self.assertTrue(callable(getattr(api,'scene_record',None)),
                        'actual flexible result already carries dt_seconds; scene merge must verify it')
        row=api.scene_record(metadata,scene_id='flex-test',model='flexible',variant='original',
            phase_index=0,dt_seconds=.0005,trajectory_sha256='test-only',actual_steps=len(result['records'])-1)
        self.assertEqual(row['dt_seconds'],.0005)
        self.assertEqual(row['actual_steps'],2)
        self.assertEqual(row['schema'],'flexible-bearing-actual-trajectory-v1')
        self.assertEqual(row['support_protocol'],result['support_protocol'])
        self.assertEqual(row['loaded_inertia_kg_m2'],result['loaded_inertia_kg_m2'])
        self.assertNotIn('id',metadata)

    def test_worker_scene_identity_conflicts_are_rejected(self):
        import flower_run_physics as api
        self.assertTrue(callable(getattr(api,'scene_record',None)))
        for field,wrong in [('dt_seconds',.001),('dt',.001),('id','different'),('actual_steps',3),
                            ('model','ideal'),('variant','optimized'),('phase_index',1),
                            ('trajectory_sha256','wrong')]:
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'SCENE_IDENTITY_MISMATCH'):
                api.scene_record({field:wrong,'diagnostic':'preserved'},scene_id='flex-test',model='flexible',
                    variant='original',phase_index=0,dt_seconds=.0005,trajectory_sha256='test-only',actual_steps=2)

    def setup_saved_records(self, root):
        # Attachment-consumption fixtures only, not real native certifications.
        folder,material,legacy=self.setup_records(root)
        loops=[[[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]]]
        problem=dict(design_region_loops=loops,extrusion_interval=[0.,1.],
                     extrusion_volume_density=1.,fixed_extrusion_layers=[])
        write(folder/'input.json',dict(meters_per_model_unit=1.,pla_density_kg_m3=1.))
        write(folder/'flower-cases.json',dict(cases=[dict(id='flower',planar_problem=problem)]))
        write(folder/'input-sha256.json',{'input.json':digest(folder/'input.json'),
                                       'flower-cases.json':digest(folder/'flower-cases.json')})
        for i in (1,2):
            path=folder/f'solve-job-0{i}/output/application.json'
            data=json.loads(path.read_text());case=data['cases'][0]
            case.update(id='flower',planar_problem=problem,solver_family='PLANAR_DOGLEG_V1',mode='FIDGET_SPINNER')
            case['attempts'][0]['solver_result']=dict(center=[1.,0.],radius=1.)
            write(path,data)
            check_path=folder/f'check-job-0{i}/output/scientific.json'
            checked=json.loads(check_path.read_text())
            checked.update(schema='saved-planar-material-moment-audit-v1',
                           input=dict(path=str(path),sha256=digest(path)))
            write(check_path,checked)
        obj=material/'retained.obj';obj.write_text('# isolated identity fixture, not simulation geometry\n')
        record=json.loads((material/'output-preparation.json').read_text())
        record.update(derived_geometry_policy='BOUNDED_SHARED_RATIONAL_ARRANGEMENT_V1',
                      input_sha256=digest(folder/'input.json'),
                      application_sha256=digest(folder/'solve-job-01/output/application.json'),
                      files={'retained':dict(file=obj.name,sha256=digest(obj))},
                      rotor_layers={},clamp_layers={},circle_center=[1.,0.],circle_radius=1.)
        write(material/'output-preparation.json',record)
        props=dict(mass_kg=1.,center_m=[0.,0.,.5],inertia_kg_m2=np.eye(3).tolist())
        proof=dict(schema='paper-flower-saved-independent-check-v1',full_material_gate='PASSED',
                   simulation_admission=True,original=props,derived_optimized=props,clamp=props,
                   selected_seed=0,circle_center=[1.,0.],circle_radius=1.,
                   partitions={name:dict(completed=True,complete_partition=True,no_positive_volume_overlap=True)
                               for name in ('editable','middle','original')},
                   bearing_subsets={name:dict(completed=True,a_subset_b=True) for name in ('outer-race','bearing-seat')},
                   material_pairs=[],moment_errors={name:[0.]*10 for name in ('original','optimized','authority')},
                   sources={str(p):digest(p) for p in [folder/'input.json',folder/'flower-cases.json',
                            folder/'input-sha256.json',folder/'solve-job-01/output/application.json',
                            material/'output-preparation.json']})
        saved=legacy/'saved-independent-check.json';write(saved,proof)
        return folder,material,legacy,saved

    def test_shared_policy_cannot_borrow_legacy_admission(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            # The legacy certificate is deliberately still PASSED.
            record=json.loads((material/'output-preparation.json').read_text())
            write(legacy/'exact-audit.json',dict(preparation_sha256=digest(material/'output-preparation.json'),
                full_material_gate='PASSED',simulation_admission=True))
            with self.assertRaisesRegex(ValueError,'SAVED_MATERIAL_CHECK_REQUIRED'):
                admit(material,legacy,folder)

    def test_saved_gate_rejects_preflight_and_missing_source_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            original=json.loads(saved.read_text())
            for changed in [dict(original,schema='paper-flower-saved-independent-preflight-v1'),
                            dict(original,sources={})]:
                write(saved,changed)
                with self.assertRaises(ValueError):admit(material,legacy,folder,saved_material_check=saved)

    def test_saved_gate_uses_reloaded_properties_and_rejects_unbalanced_material(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            record,proof=admit(material,legacy,folder,saved_material_check=saved)
            self.assertEqual(record['optimized']['center_m'],[0.,0.,.5])
            self.assertTrue(proof['independent_saved_balance']['passed'])
            changed=json.loads(saved.read_text());changed['derived_optimized']['center_m'][0]=.1
            write(saved,changed)
            with self.assertRaisesRegex(ValueError,'SAVED_MATERIAL_BALANCE'):
                admit(material,legacy,folder,saved_material_check=saved)

    def test_saved_gate_rejects_changed_obj_and_stale_preparation(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            (material/'retained.obj').write_text('# changed\n')
            with self.assertRaisesRegex(ValueError,'PHYSICS_MESH'):
                admit(material,legacy,folder,saved_material_check=saved)
            record=json.loads((material/'output-preparation.json').read_text());record['files']['retained']['sha256']=digest(material/'retained.obj')
            write(material/'output-preparation.json',record)
            with self.assertRaisesRegex(ValueError,'SAVED_CHECK_SOURCE'):
                admit(material,legacy,folder,saved_material_check=saved)

    def test_saved_gate_rejects_false_native_proof_even_with_pass_flags(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            proof=json.loads(saved.read_text())
            proof['partitions']['editable']['no_positive_volume_overlap']=False
            write(saved,proof)
            with self.assertRaisesRegex(ValueError,'SAVED_CHECK_GEOMETRY_PROOF'):
                admit(material,legacy,folder,saved_material_check=saved)

    def test_saved_gate_requires_both_independent_application_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            path=folder/'check-job-02/output/scientific.json'
            checked=json.loads(path.read_text());checked['input']['sha256']='0'*64;write(path,checked)
            with self.assertRaisesRegex(ValueError,'APPLICATION_CHECK_SOURCE'):
                admit(material,legacy,folder,saved_material_check=saved)

    def test_saved_gate_rejects_missing_overlap_proof(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            record=json.loads((material/'output-preparation.json').read_text())
            record['rotor_layers']['bottom']=1.
            record['files']['bottom']=record['files']['retained']
            write(material/'output-preparation.json',record)
            proof=json.loads(saved.read_text());proof['sources'][str(material/'output-preparation.json')]=digest(material/'output-preparation.json')
            write(saved,proof)
            with self.assertRaisesRegex(ValueError,'OVERLAP_PROOF'):
                admit(material,legacy,folder,saved_material_check=saved)

    def test_saved_gate_rejects_mismatched_circle_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,legacy,saved=self.setup_saved_records(Path(directory))
            proof=json.loads(saved.read_text());proof['circle_radius']=2.;write(saved,proof)
            with self.assertRaisesRegex(ValueError,'SAVED_CHECK_CIRCLE'):
                admit(material,legacy,folder,saved_material_check=saved)

    def setup_records(self, root):
        folder, material, audit = [root/x for x in ('case', 'material', 'audit')]
        write(folder/'input.json', {'frozen': True})
        application = {'cases': [{'scientific_hash': 'same', 'attempts': [dict(seed=0, status='CONVERGED', run=0,restart=0)]}]}
        for i in (1,2):
            write(folder/f'solve-job-0{i}/output/application.json', application)
            write(folder/f'check-job-0{i}/output/scientific.json', {'cases':[{'attempts':[dict(run=0,restart=0,moments_match=True,balance_constraints_pass=True)]}]})
        record = dict(files={}, selected_seed=0, input_sha256=digest(folder/'input.json'),
            application_sha256=digest(folder/'solve-job-01/output/application.json'))
        write(material/'output-preparation.json', record)
        write(audit/'exact-audit.json', dict(preparation_sha256=digest(material/'output-preparation.json'),
            full_material_gate='PASSED', simulation_admission=True))
        return folder, material, audit

    def test_unrelated_application_cannot_borrow_material_admission(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,audit=self.setup_records(Path(directory))
            for i in (1,2):
                path=folder/f'solve-job-0{i}/output/application.json'
                value=json.loads(path.read_text());value['different_input']=True;write(path,value)
            with self.assertRaisesRegex(ValueError,'APPLICATION_IDENTITY'):
                admit(material,audit,folder)

    def test_matching_material_and_solver_records_are_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,audit=self.setup_records(Path(directory))
            record,proof=admit(material,audit,folder)
            self.assertEqual(record['input_sha256'],digest(folder/'input.json'))
            self.assertTrue(proof['simulation_admission'])

    def test_changed_frozen_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder,material,audit=self.setup_records(Path(directory))
            write(folder/'input.json',{'frozen':False})
            with self.assertRaisesRegex(ValueError,'INPUT_IDENTITY'):admit(material,audit,folder)

    def test_changed_render_obj_is_rejected_before_visual_loading(self):
        from flower_render_physics import validate_material_meshes
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory);obj=folder/'part.obj';obj.write_text('v 0 0 0\n')
            record={'files':{'part':dict(file='part.obj',sha256=digest(obj))}}
            validate_material_meshes(folder,record)
            obj.write_text('v 1 0 0\n')
            with self.assertRaisesRegex(ValueError,'MESH'):validate_material_meshes(folder,record)

    def test_clamp_mass_cannot_be_replaced_by_legacy_geometry(self):
        from flower_run_physics import validate_clamp_mass_agreement
        masses={'cap':{'mass_kg':.001},'shaft':{'mass_kg':.002}}
        validate_clamp_mass_agreement({'clamp':{'mass_kg':.003}},masses)
        with self.assertRaisesRegex(ValueError,'CLAMP_MASS'):
            validate_clamp_mass_agreement({'clamp':{'mass_kg':.004}},masses)

    def test_timestep_comparison_detects_equal_peaks_but_different_motion(self):
        from flower_run_physics import compare_time_step_trajectories
        coarse=dict(t=np.array([0.,.5,1.]),support_xy_m=np.array([[0.,0.],[1.,0.],[0.,0.]]),
            angle=np.array([0.,1.,2.]),omega=np.array([2.,2.,2.]),energy_J=np.ones(3),
            relative_energy_balance_error=np.zeros(3),reaction=np.zeros((3,6)))
        fine={k:np.repeat(v,2,axis=0)[:5].copy() for k,v in coarse.items()}
        fine['t']=np.linspace(0,1,5);fine['support_xy_m']=np.array([[0.,0.],[0.,0.],[0.,0.],[1.,0.],[0.,0.]])
        fine['omega']=np.array([2.,2.,3.,2.,2.]);fine['energy_J']=np.array([1.,1.,1.1,1.,1.])
        result=compare_time_step_trajectories(coarse,fine,'flexible')
        self.assertGreater(result['maximum_support_position_difference_m'],0)
        self.assertGreater(result['maximum_angular_speed_difference_rad_s'],0)
        self.assertGreater(result['maximum_energy_difference_J'],0)

    def test_timestep_comparison_rejects_misaligned_samples(self):
        from flower_run_physics import compare_time_step_trajectories
        coarse={'t':np.array([0.,.5,1.])}
        fine={'t':np.array([0.,.3,.6,.9,1.2])}
        with self.assertRaisesRegex(ValueError,'COMMON_GRID'):
            compare_time_step_trajectories(coarse,fine,'flexible')


if __name__=='__main__':unittest.main()
