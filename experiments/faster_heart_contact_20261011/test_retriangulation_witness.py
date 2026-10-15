import copy
import importlib.util
from pathlib import Path
import unittest


class RetriangulationWitnessTests(unittest.TestCase):
    def check(self, a, b):
        path = Path(__file__).with_name('retriangulation_witness.py')
        self.assertTrue(path.exists(), 'exact oriented-chain witness is not implemented')
        spec = importlib.util.spec_from_file_location('witness', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.compare_retriangulation(a, b)

    def pair(self):
        a = dict(material_vertices_target_authority=[['0','0','0'], ['1','0','0'],
                 ['1','1','0'], ['0','1','0']], material_vertex_provenance=['a','b','c','d'],
                 material_faces=[[0,1,2],[0,2,3]], material_birth_faces=[0,0],
                 design_vertices_target_authority=[['-1','-1','0'],['4','-1','0'],['-1','4','0']],
                 design_faces=[[0,1,2]], source_vertices_target_authority=[], source_faces=[])
        b = copy.deepcopy(a)
        b['material_faces'] = [[0,1,3],[1,2,3]]
        return a,b

    def test_changed_diagonal_with_same_oriented_planar_region(self):
        a,b = self.pair()
        self.assertTrue(self.check(a,b)['proved'])

    def test_nonzero_exact_displacement_cannot_be_hidden_by_float_rounding(self):
        a,b = self.pair()
        b['material_vertices_target_authority'][2][2] = '1/'+str(2**1100)
        self.assertFalse(self.check(a,b)['proved'])

    def test_reversed_face_is_not_an_equivalent_oriented_boundary(self):
        a,b = self.pair()
        b['material_faces'][0].reverse()
        self.assertFalse(self.check(a,b)['proved'])

    def test_missing_face_or_ambiguous_vertex_provenance_rejected(self):
        a,b = self.pair()
        b['material_faces'].pop(); b['material_birth_faces'].pop()
        self.assertFalse(self.check(a,b)['proved'])
        a,b = self.pair()
        b['material_vertex_provenance'][1] = 'a'
        self.assertFalse(self.check(a,b)['proved'])


if __name__ == '__main__':
    unittest.main()
