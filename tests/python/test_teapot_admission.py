import copy
import tempfile
import unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/teapot_spin2014'))
import admission
from teapot import file_hashes


class AdmissionTests(unittest.TestCase):
    def test_rejects_wrong_case_or_objective(self):
        self.assertTrue(callable(getattr(admission,'check_case_context',None)), 'missing case gate')
        from make_it_stand.verification.spin_common import digest
        # Metadata fixture only; it is never used to mint a sealed artifact.
        native={'problem_digest':'input-a','solver_settings':{'objective':'TOP_INERTIA_ONLY'},
            'source_settings':{'maximum_depth':6},'boundary_settings':{},'density':1.,'minimum_outer_wall':.015,
            'status':'BINARY_MOMENT_VALIDATED','artifact_scientific_digest':'artifact-a'}
        config={'strategy':'inertia','solver':native['solver_settings'],'geometry':native['source_settings'],
            'boundary':{},'minimum_wall':.015,'meters_per_model_unit':.025,'density_kg_m3':1200.}
        summary={'strategy':'inertia','problem_digest':'input-a','scientific_digest':digest(native),
            'status':'BINARY_MOMENT_VALIDATED','sealed':True}
        admission.check_case_context(native,config,summary,'inertia','input-a')
        for strategy,input_id in [('lowering','input-a'),('inertia','input-b')]:
            with self.assertRaises(ValueError):admission.check_case_context(native,config,summary,strategy,input_id)
        tampered=copy.deepcopy(native);tampered['artifact_scientific_digest']='artifact-b'
        with self.assertRaises(ValueError):admission.check_case_context(tampered,config,summary,'inertia','input-a')

    def test_collider_must_belong_to_same_problem(self):
        self.assertTrue(callable(getattr(admission,'check_collider_context',None)), 'missing collider gate')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'mesh';p.write_bytes(b'mesh')
            manifest={'input_digest':'right','outer_hash':file_hashes([p]),'source_hashes':{},
                'profiles':{'default':{'hashes':file_hashes([p]),'near_upright_passed':True}}}
            admission.check_collider_context(manifest,'right')
            with self.assertRaises(ValueError):admission.check_collider_context(manifest,'wrong')

    def test_replaced_trajectory_rejected_by_producer_hash(self):
        self.assertTrue(callable(getattr(admission,'check_video_context',None)), 'missing video gate')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'trajectory.json';p.write_bytes(b'actual trajectory')
            report={'problem_digest':'right','source_unchanged':True,'trajectory_hashes':file_hashes([p]),
                'runs':[{'name':'trajectory','status':'RECORDED','path':str(p.resolve())}]}
            admission.check_video_context(report,'right',p)
            with self.assertRaises(ValueError):admission.check_video_context(report,'wrong',p)
            p.write_bytes(b'replaced trajectory')
            with self.assertRaises(ValueError):admission.check_video_context(report,'right',p)


if __name__=='__main__':unittest.main()
