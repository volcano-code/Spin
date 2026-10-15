"""Saved square-objective checks; not material/manufacturing certification."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
import numpy as np
from test_spin2014_quadric_interface import box
from make_it_stand.spin import (SpinItFasterRequest, Spin2014QuadricObjectiveConfig,
    solve_spin2014_static_quadric_hollowing)
from make_it_stand.verification.faster_material import _json_value
from make_it_stand.verification.static_material import audit_spin2014_static_directory


class SavedSquaredMaterial(unittest.TestCase):
    def saved(self, root, length_scale=1.):
        outer,faces=box(0.,1.)
        outer=tuple((x,y,4*z) for x,y,z in outer)
        design,inner=box(.1,.9)
        request=SpinItFasterRequest(outer,faces,design,inner,
            mode="STATIC_BALANCE",material_density=3.,target_origin=(.5,.5,0.))
        objective=Spin2014QuadricObjectiveConfig(objective="STATIC_LOWERING",inertia_weight=0.,
            solver_length_unit_in_declared_units=length_scale)
        run=solve_spin2014_static_quadric_hollowing(request,objective,fixed_material_density=2.)
        inputs=_json_value(dict(request=asdict(request),objective=asdict(objective),
            fixed_material_density=2.,actual_support_hull=[[0.,0.,0.],[1.,0.,0.],
                [1.,1.,0.],[0.,1.,0.]]))
        (root/"frozen-input.json").write_text(json.dumps(inputs))
        (root/"native-run.json").write_text(json.dumps(_json_value(asdict(run))))
        np.savez(root/"input-geometry.npz",outer_vertices=request.outer_vertices,outer_faces=request.outer_faces,
            design_vertices=request.design_vertices,design_faces=request.design_faces)
        return inputs

    def test_independent_saved_square_units_density_and_support(self):
        # Break caught: calling the native integrator, omitting fixed density,
        # accepting h instead of h², or claiming optimization from zero updates.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.saved(root,3.)
            audit=audit_spin2014_static_directory(root)
            from make_it_stand.verification.run_evidence import _nonempty_object
            _nonempty_object(audit,"independent static report")
            h=16.256/8.512
            self.assertAlmostEqual(audit["independent_total_raw_moments"][0],8.512,places=12)
            self.assertAlmostEqual(audit["independent_objective"],9*h*h,places=11)
            self.assertTrue(audit["terminal_numerics_qualified"])
            self.assertTrue(audit["support"]["inside_support"])
            self.assertFalse(audit["nontrivial_optimization"])
            self.assertFalse(audit["material_partition_verified"])
            self.assertFalse(audit["manufacturing_approved"])

    def test_changed_density_or_squared_formula_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inputs=self.saved(root)
            record=json.loads((root/"native-run.json").read_text())
            record["squared_height"]["objective"]=16.256/8.512
            (root/"native-run.json").write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError,"objective"):
                audit_spin2014_static_directory(root)
            self.saved(root)
            inputs["fixed_material_density"]=4.
            (root/"frozen-input.json").write_text(json.dumps(inputs))
            with self.assertRaisesRegex(ValueError,"density|moment"):
                audit_spin2014_static_directory(root)

    def test_changed_geometry_or_forged_support_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inputs=self.saved(root)
            inputs["request"]["outer_vertices"][0][0]=.01
            (root/"frozen-input.json").write_text(json.dumps(inputs))
            with self.assertRaisesRegex(ValueError,"geometry"):
                audit_spin2014_static_directory(root)
            inputs=self.saved(root)
            inputs["actual_support_hull"]=[[0.,0.,.1],[1.,0.,.1],[1.,1.,.1],[0.,1.,.1]]
            (root/"frozen-input.json").write_text(json.dumps(inputs))
            with self.assertRaisesRegex(ValueError,"support"):
                audit_spin2014_static_directory(root)

    def test_existing_checker_cli_records_independent_process_without_native(self):
        # Break caught: a new static result being checked by the old h-objective
        # parser or importing the optimizer to compute its own expectations.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.saved(root)
            repository=Path(__file__).resolve().parents[2]
            output=root/"check"
            process=subprocess.run([sys.executable,str(repository/
                "experiments/spin_it_reproduction/faster/check_saved_static_material.py"),
                "--spin2014-function-directory",str(root),"--output",str(output)],
                capture_output=True,text=True,timeout=60)
            self.assertEqual(process.returncode,0,process.stdout+process.stderr)
            report=json.loads((output/"scientific.json").read_text())
            self.assertEqual(report["schema"],"saved-spin2014-static-square-audit-v1")
            self.assertTrue(report["terminal_numerics_qualified"])
            self.assertFalse(report["material_partition_verified"])
            message=json.loads(process.stdout.strip())
            self.assertFalse(message["native_solver_loaded"])
            self.assertTrue((output/"evidence.json").is_file())


if __name__=="__main__":unittest.main()
