"""Exercise the native settings/trace boundary, not hand-written pass flags."""
from dataclasses import asdict, replace
import json
import unittest
from test_source_spin_material_bridge import problem, settings
from make_it_stand.verification import source_spin_material as bridge
import _make_it_stand as native


class PaperTerminationContract(unittest.TestCase):
    def test_explicit_native_policy_and_kkt_scope(self):
        p=problem();cfg=settings(bridge)
        args=dict(source_settings=asdict(cfg['source_settings']),
                  solver_settings=asdict(cfg['solver_settings']),
                  boundary_settings=asdict(cfg['boundary_settings']),
                  density=1.,minimum_outer_wall=.5)
        args['solver_settings']['termination_policy']='PAPER_POST_TRANSITION_V1'
        try:
            result=native.solve_source_spin2014_material(
                *p.original.arrays(),*p.design_regions[0].mesh.arrays(),**args)
        except ValueError as exc:
            self.fail('native must support the explicit paper termination policy: '+str(exc))
        facts=result.diagnostics;solve=facts['solve']
        self.assertEqual(solve['termination_policy'],'PAPER_POST_TRANSITION_V1')
        self.assertEqual(
            facts['settings']['solver_settings'],
            json.loads(json.dumps(args['solver_settings'])),
        )
        self.assertEqual(solve['status'],'BINARY_MOMENT_VALIDATED')
        self.assertEqual(len(solve['passes']),1)
        self.assertFalse(solve['terminal_relaxation_kkt_applies'])
        last=solve['passes'][-1]
        self.assertTrue(last['final_kkt']['valid'])
        self.assertEqual(last['kkt_partition_hash'],last['before_hash'])
        self.assertNotEqual(last['kkt_partition_hash'],last['after_hash'])
        self.assertTrue(last['paper_resolution_complete'])
        self.assertEqual(last['output_fractional_count'],0)
        # The high-level bridge must reintegrate and seal the actual material,
        # not reject a valid binary object merely because its relaxed inherited
        # partition lacks a new KKT, nor use that KKT as a substitute for a seal.
        cfg['solver_settings']=replace(cfg['solver_settings'],termination_policy='PAPER_POST_TRANSITION_V1')
        sealed=bridge.solve_and_seal_source_spin2014(p,**cfg)
        self.assertIsNotNone(sealed.artifact)
        self.assertIsNone(sealed.seal_error)
        self.assertEqual(sealed.identity['native_run']['solve'],solve)
        self.assertAlmostEqual(sealed.artifact.facts['physical']['mass'],503/64)
        args['solver_settings']['termination_policy']='FORGED_PASS'
        with self.assertRaises(ValueError):
            native.solve_source_spin2014_material(
                *p.original.arrays(),*p.design_regions[0].mesh.arrays(),**args)

    def test_legacy_omission_is_not_a_paper_upgrade(self):
        p=problem();cfg=settings(bridge)
        args=dict(source_settings=asdict(cfg['source_settings']),
                  solver_settings=asdict(cfg['solver_settings']),
                  boundary_settings=asdict(cfg['boundary_settings']),
                  density=1.,minimum_outer_wall=.5)
        explicit=native.solve_source_spin2014_material(
            *p.original.arrays(),*p.design_regions[0].mesh.arrays(),**args).diagnostics
        args['solver_settings'].pop('termination_policy')
        legacy=native.solve_source_spin2014_material(
            *p.original.arrays(),*p.design_regions[0].mesh.arrays(),**args).diagnostics
        self.assertEqual(explicit,legacy)
        self.assertEqual(legacy['solve']['termination_policy'],'UNCHANGED_PARTITION_KKT_V1')
        self.assertTrue(legacy['solve']['terminal_relaxation_kkt_applies'])
        for bad in ('FORGED_PASS',True,None,[]):
            with self.assertRaises(ValueError):
                bridge.SourceSolverSettings(termination_policy=bad)


if __name__=='__main__':unittest.main()
