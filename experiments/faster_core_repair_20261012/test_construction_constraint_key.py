import unittest
import cherchi_core_map as core


class ConstraintKeyTests(unittest.TestCase):
    def key(self,kind,indices):
        fn=getattr(core,'construction_constraint_key',None)
        self.assertTrue(callable(fn),'typed construction constraint key missing')
        return fn(dict(point_type=kind,operands=[dict(status='UNIQUE',source_vertices=[
            dict(input_object=0,input_vertex=i)]) for i in indices]))

    def test_lpi_line_direction_and_plane_corner_order_do_not_change_constraints(self):
        self.assertEqual(self.key('LPI',[0,1,2,3,4]),self.key('LPI',[1,0,4,2,3]))

    def test_changed_lpi_plane_changes_constraints(self):
        self.assertNotEqual(self.key('LPI',[0,1,2,3,4]),self.key('LPI',[0,1,2,3,5]))

    def test_tpi_plane_order_does_not_change_constraints(self):
        self.assertEqual(self.key('TPI',list(range(9))),self.key('TPI',[6,7,8,3,4,5,0,1,2]))

    def replay(self,vertices):
        fn=getattr(core,'replay_lpi_constraint',None)
        self.assertTrue(callable(fn),'independent exact LPI replay missing')
        record=dict(point_type='LPI',operands=[dict(status='UNIQUE',source_vertices=[
            dict(input_object=0,input_vertex=i)]) for i in range(5)])
        return fn(record,[vertices])

    def test_independent_lpi_hand_plane_intersection(self):
        self.assertEqual(self.replay([(0,0,-1),(0,0,1),(0,0,0),(1,0,0),(0,1,0)]),(0,0,0))

    def test_parallel_lpi_cannot_fabricate_point(self):
        with self.assertRaisesRegex(ValueError,'SINGULAR_LPI'):
            self.replay([(0,0,1),(1,0,1),(0,0,0),(1,0,0),(0,1,0)])
