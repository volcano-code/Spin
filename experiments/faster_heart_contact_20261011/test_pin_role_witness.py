import tempfile
from pathlib import Path
import unittest


class PinRoleWitnessTests(unittest.TestCase):
    def module(self):
        import importlib.util
        self.assertIsNotNone(importlib.util.find_spec('pin_role_witness'),
                             'pin source parsing and exact difference witness entry must exist')
        import pin_role_witness
        self.assertTrue(callable(getattr(pin_role_witness,'load_pin_boundary',None)))
        return pin_role_witness

    def test_obj_binary64_coordinates_are_preserved_without_decimal_authority(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pin.obj'
            path.write_text('v 0.1 0 0\nv 0 1 0\nv 0 0 1\nf 1 2 3\n')
            packet=self.module().load_pin_boundary(path)
            self.assertEqual(packet['vertices'][0][0],'3602879701896397/36028797018963968')
            self.assertEqual(packet['faces'],[[0,1,2]])

    def test_obj_nonfinite_negative_or_nontriangle_faces_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pin.obj'
            for data in ('v nan 0 0\nf 1 1 1\n','v 0 0 0\nf -1 1 1\n',
                         'v 0 0 0\nf 1 1 1 1\n'):
                path.write_text(data)
                with self.assertRaises(ValueError):self.module().load_pin_boundary(path)

    def test_matching_vertices_do_not_prove_entity_containment(self):
        packet=dict(vertices=[['0/1','0/1','0/1'],['1/1','0/1','0/1']],faces=[])
        report=self.module().coordinate_correspondence(packet,packet)
        self.assertEqual(report['source_pin_vertices_found_exactly'],2)
        self.assertFalse(report['containment_proved'])

    def test_entire_corresponding_triangle_plane_departure_is_exact_not_sampling(self):
        source=dict(vertices=[['0/1','0/1','0/1'],['1/1','0/1','0/1'],['0/1','1/1','0/1']],faces=[[0,1,2]])
        fixed=dict(vertices=source['vertices'],faces=[[0,1,2]])
        derived=dict(vertices=[['0/1','0/1','-1/3'],*source['vertices'][1:]],faces=[[0,1,2]])
        report=self.module().corresponding_pin_planes(source,fixed,derived)
        self.assertEqual(report['matched_authority_fixed_faces'],1)
        self.assertEqual(report['departed_faces'][0]['signed_plane_vertex_values_exact'],['-1/3','0/1','0/1'])
        self.assertFalse(report['containment_proved'])

if __name__=='__main__':unittest.main()
