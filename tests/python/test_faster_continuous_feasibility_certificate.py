"""The frozen offset-axis fixture has an exact continuous moment-capacity bound."""
from fractions import Fraction
from dataclasses import replace
import importlib
import json
from pathlib import Path
import tempfile
import unittest


class FasterContinuousFeasibilityCertificateTests(unittest.TestCase):
    def setUp(self):
        self.repo = Path(__file__).resolve().parents[2]
        self.fixture = self.repo / "experiments/spin_it_reproduction/faster/fixtures/cherchi-rotated-box-fixed-shell-v2/problem.json"
        self.checker = importlib.import_module("make_it_stand.verification.faster_infeasibility")

    def _nested_box_manifest(self, root, *, sign=1, target=-1, density=1):
        # Literal oriented box topology; expectations below come from volumes
        # and centroids, not the production integrator or capacity checker.
        faces = ((0,2,1),(0,3,2),(4,5,6),(4,6,7),(0,1,5),(0,5,4),
                 (1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,0,4),(3,4,7))
        def box(lo, hi):
            vertices = ((lo[0],lo[1],lo[2]),(hi[0],lo[1],lo[2]),
                (hi[0],hi[1],lo[2]),(lo[0],hi[1],lo[2]),
                (lo[0],lo[1],hi[2]),(hi[0],lo[1],hi[2]),
                (hi[0],hi[1],hi[2]),(lo[0],hi[1],hi[2]))
            return tuple((sign*x,y,z) for x,y,z in vertices)
        oriented = faces if sign==1 else tuple((a,c,b) for a,b,c in faces)
        outer = box((-4,-2,-2),(2,2,2)); inner = box((0,-1,-1),(1,1,1))
        shell_faces = oriented + tuple((a+8,c+8,b+8) for a,b,c in oriented)
        for name, vertices, triangles in (("outer",outer,oriented),("design",inner,oriented),
                                          ("shell",outer+inner,shell_faces)):
            lines = ["v " + " ".join(map(str,v)) for v in vertices]
            lines += ["f " + " ".join(str(i+1) for i in f) for f in triangles]
            (root/f"{name}.obj").write_text("\n".join(lines)+"\n")
        data = {"schema":"faster-cherchi-problem-manifest-v1", "outer_mesh":"outer.obj",
            "fixed_regions":[{"mesh":"shell.obj","density":1}],
            "design_regions":[{"mesh":"design.obj","density":density}],
            "target":{"mode":"YO_YO","axis":[sign,0,0],"origin":[sign*target,0,0]}}
        manifest = root/"problem.json"; manifest.write_text(json.dumps(data))
        return manifest

    def test_exterior_target_capacity_cannot_reject_full_material_witness(self):
        # V_outer=96, centroid=(-1,0,0); V_design=4, centroid=(1/2,0,0).
        # Relative shell moment=-6. Retaining the complete design contributes
        # +6, so zero is attainable; old width-squared formula gave only +2.
        for sign in (1,-1):
            with self.subTest(sign=sign), tempfile.TemporaryDirectory() as temp:
                report = self.checker.verify_axis_aligned_yoyo_moment_capacity(
                    self._nested_box_manifest(Path(temp), sign=sign))
                self.assertEqual(report.status, "UNKNOWN")
                cert = report.continuous_certificate
                self.assertEqual(cert.fixed_relative_first_moment, Fraction(-6))
                self.assertEqual(cert.maximum_editable_positive_moment, Fraction(6))
                self.assertEqual(cert.total_moment_upper_bound, Fraction(0))

    def test_independent_certificate_rejects_old_underestimated_capacity(self):
        from make_it_stand.verification.feasibility import ContinuousMomentCapacityCertificate
        cert = ContinuousMomentCapacityCertificate(
            constrained_axis=0,direction_sign=1,target_coordinate=Fraction(-1),
            fixed_density=Fraction(1),fixed_volume=Fraction(92),fixed_first_moment=Fraction(-98),
            fixed_relative_first_moment=Fraction(-6),editable_density=Fraction(1),
            box_lower=(Fraction(0),Fraction(-1),Fraction(-1)),
            box_upper=(Fraction(1),Fraction(1),Fraction(1)),
            maximum_editable_positive_moment=Fraction(2),total_moment_upper_bound=Fraction(-4),
            input_sha256="a"*64,asset_sha256=(("mesh", "b"*64),))
        self.assertFalse(cert.valid, "old width-only capacity is not an upper bound")
        # At half density the correct bound is -6+3=-3, genuinely impossible.
        self.assertTrue(replace(cert,editable_density=Fraction(1,2),
            maximum_editable_positive_moment=Fraction(3),total_moment_upper_bound=Fraction(-3)).valid)

    def test_capacity_has_exact_outside_inside_and_boundary_values(self):
        # Integral over [0,1] times transverse area4. All expectations literal.
        for target, capacity in ((-1,6),(0,2),(0.5,Fraction(1,2)),(1,0),(2,0)):
            for sign in (1,-1):
                with self.subTest(target=target,sign=sign), tempfile.TemporaryDirectory() as temp:
                    report = self.checker.verify_axis_aligned_yoyo_moment_capacity(
                        self._nested_box_manifest(Path(temp),sign=sign,target=target))
                    self.assertEqual(report.continuous_certificate.maximum_editable_positive_moment,
                                     Fraction(capacity))

    def test_exterior_target_still_certifies_insufficient_material_density(self):
        for sign in (1,-1):
            with self.subTest(sign=sign), tempfile.TemporaryDirectory() as temp:
                report = self.checker.verify_axis_aligned_yoyo_moment_capacity(
                    self._nested_box_manifest(Path(temp),sign=sign,density=0.5))
                self.assertEqual(report.status,"INFEASIBLE_CERTIFIED")
                self.assertEqual(report.continuous_certificate.maximum_editable_positive_moment,Fraction(3))
                self.assertEqual(report.continuous_certificate.total_moment_upper_bound,Fraction(-3))

    def test_offset_fixed_shell_is_certified_infeasible_from_meshes(self):
        report = self.checker.verify_axis_aligned_yoyo_moment_capacity(self.fixture)

        self.assertEqual(report.status, "INFEASIBLE_CERTIFIED")
        self.assertEqual(report.scope, "CONTINUOUS_GEOMETRIC_MOMENT_CAPACITY")
        certificate = report.continuous_certificate
        self.assertTrue(certificate.valid)
        self.assertEqual(certificate.fixed_relative_first_moment, Fraction(-21, 2))
        self.assertEqual(certificate.maximum_editable_positive_moment, Fraction(9, 4))
        self.assertEqual(certificate.total_moment_upper_bound, Fraction(-33, 4))
        self.assertEqual(certificate.constrained_axis, 0)
        self.assertEqual(certificate.input_sha256, report.input_hash)

    def test_centered_target_does_not_reuse_offset_certificate(self):
        manifest = json.loads(self.fixture.read_text())
        manifest["target"]["origin"] = [0.0, 0.0, 0.0]
        with tempfile.TemporaryDirectory() as temp:
            temp_root = Path(temp)
            for name in ("outer.obj", "design.obj", "shell.obj"):
                (temp_root / name).write_bytes((self.fixture.parent.parent / "cherchi-rotated-box-fixed-shell-v1" / name).read_bytes())
            manifest["outer_mesh"] = "outer.obj"
            manifest["design_regions"][0]["mesh"] = "design.obj"
            manifest["fixed_regions"][0]["mesh"] = "shell.obj"
            changed = temp_root / "problem.json"
            changed.write_text(json.dumps(manifest))
            report = self.checker.verify_axis_aligned_yoyo_moment_capacity(changed)

        self.assertEqual(report.status, "UNKNOWN")
        self.assertFalse(report.continuous_certificate.valid)
        self.assertNotEqual(report.input_hash, self.checker.verify_axis_aligned_yoyo_moment_capacity(self.fixture).input_hash)

    def test_expected_input_identity_mismatch_cannot_reuse_certificate(self):
        report = self.checker.verify_axis_aligned_yoyo_moment_capacity(self.fixture, expected_input_hash="0" * 64)
        self.assertEqual(report.status, "UNKNOWN")
        self.assertIn("INPUT_IDENTITY_MISMATCH", report.reason)

    def test_centered_fixture_has_an_independent_non_degenerate_yoyo_witness(self):
        centered = self.fixture.parent.parent / "cherchi-centered-axis-fixed-shell-v1/problem.json"
        report = self.checker.verify_centered_axis_ellipsoid_witness(centered)

        self.assertEqual(report.status, "FEASIBLE_WITNESS")
        witness = report.continuous_witness
        self.assertTrue(witness.valid)
        self.assertEqual(witness.ellipsoid_semiaxes, (Fraction(1, 2), Fraction(3, 4), Fraction(3, 4)))
        self.assertEqual(witness.inertia_gap_per_mass, Fraction(1, 16))
        self.assertEqual(witness.input_sha256, report.input_hash)

    def test_offset_target_cannot_use_centered_symmetry_witness(self):
        report = self.checker.verify_centered_axis_ellipsoid_witness(self.fixture)
        self.assertEqual(report.status, "UNKNOWN")
        self.assertFalse(report.continuous_witness.valid)

    def test_certificate_is_bound_to_report_input_and_serializes_scope(self):
        report = self.checker.verify_axis_aligned_yoyo_moment_capacity(self.fixture)
        encoded = report.to_dict()
        self.assertEqual(encoded["schema"], "continuous-geometric-feasibility-v1")
        self.assertTrue(encoded["continuous_certificate"]["valid"])
        self.assertEqual(encoded["continuous_certificate"]["input_sha256"], report.input_hash)
        stale = replace(report, input_hash="f" * 64)
        self.assertEqual(stale.status, "UNKNOWN")


if __name__ == "__main__":
    unittest.main()
