"""Figure 11 volume semantics, independent of uniform density and legacy mass."""
import copy
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
import numpy as np
from test_faster_static_table_dispatch import method_case
from test_spin_it_faster_table1_runner import load_runner
from make_it_stand.spin import solve_spin_it_faster_static_plane_family

VOLUME="UNIFORM_DENSITY_VOLUME_V1"
MASS="DENSITY_WEIGHTED_MASS_V1"

class StaticVolumeTests(unittest.TestCase):
    def test_external_evaluator_rejects_wrong_fixed_density_before_geometry(self):
        from make_it_stand.spin import solve_spin_it_faster_with_external_moment_evaluator
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            original=runner.build_request(method_case(root,"PARETO_TRADEOFF",.002),asset_root=root)
            request=replace(original,static_objective_measure=VOLUME,
                fixed_shell_moments=tuple(2*x for x in original.fixed_shell_moments),
                settings=replace(original.settings,equation_profile="DUAL_MONOMIAL_RAW_V1"))
            def unexpected_geometry(_):
                raise AssertionError("unvalidated volume request reached geometry")
            with self.assertRaisesRegex(ValueError,"uniform|fixed"):
                solve_spin_it_faster_with_external_moment_evaluator(request,unexpected_geometry,
                    initial_quadrics=[(1,0,0,0,-1,-1,-1,0,0,0)])

    def test_paper_volume_solve_is_density_invariant_but_legacy_mass_is_not(self):
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); results=[]
            for density in (1.,1.7):
                case=method_case(root,"PARETO_TRADEOFF",.002)
                case["material_density"]=density
                request=replace(runner.build_request(case,asset_root=root),static_objective_measure=VOLUME)
                result=solve_spin_it_faster_static_plane_family(request,initial_coefficients=(0.,.08,-.06,-1.))
                self.assertEqual(result["status"],"CONVERGED")
                self.assertGreater(result["accepted_steps"],0)
                results.append(result)
            np.testing.assert_allclose(results[0]["plane_coefficients"],results[1]["plane_coefficients"],atol=1e-10,rtol=0)
            np.testing.assert_allclose(results[0]["total_raw_moments"],results[1]["total_raw_moments"]/1.7,atol=1e-9,rtol=0)
            legacy=solve_spin_it_faster_static_plane_family(replace(request,static_objective_measure=MASS),
                initial_coefficients=(0.,.08,-.06,-1.))
            self.assertEqual(legacy["status"],"CONVERGED")
            self.assertGreater(np.linalg.norm(legacy["plane_coefficients"]-results[1]["plane_coefficients"]),.005)

    def test_volume_rejects_unverified_nonuniform_fixed_mass(self):
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            request=replace(runner.build_request(method_case(root,"PARETO_TRADEOFF",.002),asset_root=root),
                            static_objective_measure=VOLUME)
            with self.assertRaisesRegex(ValueError,"uniform|fixed"):
                solve_spin_it_faster_static_plane_family(replace(request,
                    fixed_shell_moments=tuple(2*x for x in request.fixed_shell_moments)),
                    initial_coefficients=(0.,.08,-.06,-1.))

    def test_formal_static_pareto_uses_volume_and_checks_saved_geometry(self):
        from test_spin_it_faster_paper_figures import load_runner as load_figures
        from make_it_stand.verification.static_material import audit_static_case_material
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result=load_figures().run_pareto_sweep(method_case(root),[0.,.002,1.],asset_root=root)
            self.assertEqual(result["static_objective_measure"],VOLUME)
            for row in result["rows"]:
                self.assertEqual(row["status"],"CONVERGED")
                application=row["application"]
                self.assertEqual(application["static_problem"]["static_objective_measure"],VOLUME)
                self.assertAlmostEqual(row["material_volume"],row["material_mass"]/1.7,places=11)
                self.assertAlmostEqual(row["first_objective"],-row["material_volume"],places=11)
                self.assertTrue(audit_static_case_material(application)["attempts"][-1]["terminal_numerics_qualified"])
            # Relabeling the same terminal geometry as the legacy objective must
            # not produce a false numerical/paper pass.
            forged=copy.deepcopy(result["rows"][1]["application"])
            forged["attempts"][-1]["static_objective_measure"]=MASS
            with self.assertRaisesRegex(ValueError,"measure|input"):
                audit_static_case_material(forged)

if __name__=="__main__":unittest.main()
