from fractions import Fraction as F
import unittest
import birth_face_recipes as recipes


class BirthFaceRecipes(unittest.TestCase):
    def setUp(self):
        self.operands = [[(F(0), F(0), F(0)), (F(2), F(0), F(0)), (F(0), F(2), F(0))],
                         [(F(1), F(-1), F(-1)), (F(1), F(3), F(-1)), (F(1), F(-1), F(3))]]
        self.faces = [[(0, 1, 2)], [(0, 1, 2)]]

    def test_edge_plane_recipe_replays_motion_not_nearest_coordinate(self):
        recipe = recipes.recover_vertex_recipe((F(1), F(0), F(0)), [0, 1], self.operands, self.faces)
        moved = [self.operands[0], [(F(3, 2), y, z) for _, y, z in self.operands[1]]]
        self.assertEqual(recipes.replay_vertex_recipe(recipe, moved), (F(3, 2), F(0), F(0)))

    def test_explicit_input_vertex_is_retained(self):
        recipe = recipes.recover_vertex_recipe(self.operands[0][0], [0], self.operands, self.faces)
        self.assertEqual(recipe, ('EXPLICIT', 0, 0))

    def test_wrong_birth_face_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'BIRTH_FACE'):
            recipes.recover_vertex_recipe((F(1), F(0), F(1, 10**60)), [0, 1], self.operands, self.faces)

    def test_rank_deficient_plane_intersection_is_not_a_recipe(self):
        operands = [self.operands[0], self.operands[0]]
        with self.assertRaisesRegex(ValueError, 'UNRESOLVED'):
            recipes.recover_vertex_recipe((F(1, 2), F(1, 2), F(0)), [0, 1], operands, self.faces)


if __name__ == '__main__': unittest.main()
