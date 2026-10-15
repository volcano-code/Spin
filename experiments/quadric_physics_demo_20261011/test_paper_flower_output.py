import unittest
import numpy as np
import trimesh
from paper_flower_output import raw_mesh_moments,CIRCLE_SEGMENTS,profile_geometry


class IndependentMeshMomentsTests(unittest.TestCase):
    def test_bearing_seat_uses_every_shared_bore_node(self):
        from paper_flower_output import shared_bearing_profiles
        # Rounded cut intersection is deliberately off the old edge. Both
        # mating regions must use the identical new vertex, not weld a copy.
        bore=[[0.,0.],[0.,2.],[1.,2.],[1.0000000000000002,1.],[1.,0.]]
        seat_loops=[[[0.,0.],[1.,0.],[1.,2.],[0.,2.]],
                    [[.2,.2],[.2,1.8],[.8,1.8],[.8,.2]]]
        hole,seat=shared_bearing_profiles([[],bore],1,seat_loops)
        self.assertEqual(set(hole.exterior.coords[:-1]),set(map(tuple,bore)))
        self.assertEqual(set(seat.exterior.coords[:-1]),set(map(tuple,bore)))
        self.assertEqual(set(seat.interiors[0].coords[:-1]),set(map(tuple,seat_loops[1])))

    def test_translated_asymmetric_box_non_diagonal_second_moments(self):
        mesh=trimesh.creation.box(extents=[2.,4.,6.]);mesh.apply_translation([1.,2.,3.])
        raw=raw_mesh_moments(mesh,2.)
        expected=2*np.array([48.,48.,96.,144.,64.,256.,576.,288.,144.,96.])
        np.testing.assert_allclose(raw,expected,rtol=1e-14,atol=1e-12)

    def test_registered_arc_resolution(self):
        self.assertEqual(CIRCLE_SEGMENTS,4096)

    def test_nested_island_is_not_subtracted_with_parent_hole(self):
        loops=[[[0,0],[4,0],[4,4],[0,4]],
               [[1,1],[1,3],[3,3],[3,1]],
               [[1.5,1.5],[2.5,1.5],[2.5,2.5],[1.5,2.5]]]
        geometry=profile_geometry(loops)
        self.assertTrue(geometry.is_valid)
        self.assertEqual(geometry.area,13.)

    def test_obj_writer_preserves_near_collinear_exact_binary64_area(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        import paper_flower_output as output
        # Actual original-flower face 887: the old fixed 17-decimal-place
        # exporter changed nonzero binary64 area into exact collinearity.
        points=[[-.001972812725272114,-.004780060504351415,-.0016],
                [-.0021405100601340397,-.004943493408309201,-.0016],
                [-.0023082073949959657,-.005106926312266985,-.0016]]
        mesh=trimesh.Trimesh(vertices=points,faces=[[0,1,2]],process=False)
        self.assertEqual(output.exact_zero_area_faces(mesh),[])
        with TemporaryDirectory() as directory:
            legacy=Path(directory)/'legacy.obj'
            mesh.export(legacy,digits=17)
            legacy_mesh=trimesh.load(legacy,force='mesh',process=False)
            self.assertEqual(len(output.exact_zero_area_faces(legacy_mesh)),1)
            path=Path(directory)/'triangle.obj'
            output.write_roundtrip_obj(mesh,path)
            loaded=trimesh.load(path,force='mesh',process=False)
            self.assertEqual(output.exact_zero_area_faces(loaded),[])
            np.testing.assert_array_equal(loaded.vertices,mesh.vertices)

    def test_48_hole_cap_retriangulation_keeps_every_boundary_vertex(self):
        from pathlib import Path
        import json
        import paper_flower_output as output
        record=json.loads((Path(__file__).resolve().parent/'paper-flower-method-v1/input.json').read_text())
        polygon=profile_geometry(record['retained_original_loops'])
        mesh=output.constrained_boundary_prism(polygon,[-.08,.08],.02)
        expected={tuple(p) for loop in record['retained_original_loops'] for p in loop}
        self.assertEqual({tuple(p[:2]) for p in mesh.vertices},expected)
        self.assertTrue(mesh.is_watertight)
        mesh.apply_scale(.02)
        self.assertEqual(output.exact_zero_area_faces(mesh),[])

    def test_audit_exception_preserves_failure_attachment_and_event_prefix(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        from types import SimpleNamespace
        from unittest.mock import patch
        import json
        import paper_flower_output as output
        def reject(*args):raise ValueError('CONTROLLED_EXACT_OPERAND_FAILURE')
        native=SimpleNamespace(make_exact_material_region=reject,__file__=__file__)
        with TemporaryDirectory() as directory:
            out=Path(directory)/'derived';out.mkdir();report=Path(directory)/'audit'
            (out/'output-preparation.json').write_text(json.dumps(dict(files={},rotor_layers={})))
            output.write_roundtrip_obj(trimesh.creation.box(),out/'design.obj')
            with patch.dict('sys.modules',{'_make_it_stand':native}):
                with self.assertRaisesRegex(ValueError,'CONTROLLED_EXACT_OPERAND_FAILURE'):
                    output.exact_audit(out,report)
            failure=json.loads((report/'audit-failure.json').read_text())
            self.assertFalse(failure['simulation_admission'])
            self.assertEqual(failure['last_completed_or_started_event']['name'],'design')
            events=[json.loads(line) for line in (report/'audit-events.jsonl').read_text().splitlines()]
            self.assertEqual([e['stage'] for e in events],['audit_started','authenticate_operand','audit_failed'])

    def test_legacy_exact_audit_never_admits_from_hand_moment_flags(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        from types import SimpleNamespace
        from unittest.mock import patch
        import json,hashlib
        import paper_flower_output as output
        # The geometry dependency returns completed checks; the admission
        # decision must nevertheless require the independent saved checker.
        native=SimpleNamespace(__file__=__file__,make_exact_material_region=lambda *args:None,
            check_exact_material_partition=lambda *args:dict(completed=True,complete_partition=True,no_positive_volume_overlap=True),
            check_exact_material_subset=lambda *args:dict(completed=True,a_subset_b=True),
            relate_exact_material_regions=lambda *args:dict(completed=True,positive_volume_overlap=False))
        with TemporaryDirectory() as directory:
            out=Path(directory)/'derived';out.mkdir();report=Path(directory)/'audit'
            files={}
            for name in ['design','retained','cavity','original-middle','original-flower','bearing-bore-region','fixed-middle-boundary-rim','bearing-seat','outer-race','bottom','top']:
                path=out/(name+'.obj');output.write_roundtrip_obj(trimesh.creation.box(),path)
                files[name]=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            prep=dict(files=files,rotor_layers={'fixed-middle-boundary-rim':1,'bearing-seat':1,'outer-race':1,'bottom':1,'top':1},
                      moment_gate_pass=True,authority_independent_moment_error=[0]*10)
            (out/'output-preparation.json').write_text(json.dumps(prep))
            with patch.dict('sys.modules',{'_make_it_stand':native}):output.exact_audit(out,report)
            result=json.loads((report/'exact-audit.json').read_text())
            self.assertFalse(result['simulation_admission'])
            self.assertEqual(result['full_material_gate'],'UNKNOWN')
            self.assertTrue(result['geometry_partition_checks_passed'])


if __name__=='__main__':unittest.main()
