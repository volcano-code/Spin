"""Conservative, explicit topology retention; never upgrades a paper solve."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/heart_stem_contact'))
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/reconstructed_heart_top'))


class TopologyRetentionTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('topology_retention'),'missing audited topology repair')
        import topology_retention
        return topology_retention

    def test_edge_contact_retains_material_instead_of_adding_carving(self):
        m=self.module(); original={(0,0,0),(1,1,0),(5,5,5)}
        result,trace=m.retain_critical_cells(original)
        self.assertEqual(result,{(5,5,5)})
        self.assertEqual(set(map(tuple,trace[0]['retained_cells'])),{(0,0,0),(1,1,0)})
        self.assertEqual(original,{(0,0,0),(1,1,0),(5,5,5)})

    def test_vertex_contact_is_not_missed_by_edge_only_detection(self):
        m=self.module()
        result,_=m.retain_critical_cells({(0,0,0),(1,1,1),(5,5,5)})
        self.assertEqual(result,{(5,5,5)})

    def test_face_connected_shape_and_disjoint_cavities_remain_unchanged(self):
        m=self.module(); original={(0,0,0),(1,0,0),(1,1,0),(5,5,5)}
        result,trace=m.retain_critical_cells(original)
        self.assertEqual(result,original); self.assertEqual(trace,[])

    def test_no_cavity_is_failure_not_successful_optimization(self):
        m=self.module()
        with self.assertRaisesRegex(ValueError,'EMPTY_CAVITY'):
            m.retain_critical_cells({(0,0,0),(1,1,0)})

    def test_last_allowed_repair_pass_can_complete(self):
        m=self.module()
        result,trace=m.retain_critical_cells({(0,0,0),(1,1,0),(5,5,5)},maximum_passes=1)
        self.assertEqual(result,{(5,5,5)}); self.assertEqual(len(trace),1)

    def test_archived_solver_settings_cannot_be_relabelled(self):
        import seal_retention as s
        self.assertTrue(hasattr(s,'check_archived_configuration'),'missing configuration provenance gate')
        config=dict(geometry={'maximum_depth':6},solver={'objective':'TOP_MASS_AND_LOWERING'},
                    boundary={'maximum_input_cuboids':4000},minimum_wall=.015,
                    meters_per_model_unit=.025,density_kg_m3=1200.)
        prior=dict(source_settings=config['geometry'],solver_settings=config['solver'],
                   boundary_settings=config['boundary'],minimum_outer_wall=.015,density=1.)
        s.check_archived_configuration(config,prior)
        for field,value in [('geometry',{'maximum_depth':9}),('solver',{'objective':'STATIC_LOWERING'}),
                            ('minimum_wall',0.),('meters_per_model_unit',1.),('density_kg_m3',1.)]:
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'CONFIGURATION_MISMATCH'):
                s.check_archived_configuration({**config,field:value},prior)

    def test_exact_grid_expansion_and_resource_rejection(self):
        m=self.module()
        boxes=[dict(cell=dict(level=1,coordinates=[0,0,0]),lower=[0.,0.,0.],upper=[1.,1.,1.])]
        cells=m.expand_cuboids(boxes,depth=2,origin=(0.,0.,0.),side=2.,maximum_cells=8)
        self.assertEqual(len(cells),8)
        with self.assertRaisesRegex(ValueError,'RESOURCE_LIMIT'):
            m.expand_cuboids(boxes,depth=2,origin=(0.,0.,0.),side=2.,maximum_cells=7)
        boxes[0]['upper'][0]=1.001
        with self.assertRaisesRegex(ValueError,'GRID_GEOMETRY_MISMATCH'):
            m.expand_cuboids(boxes,depth=2,origin=(0.,0.,0.),side=2.)

    def test_surface_and_independent_box_moments_include_cross_terms(self):
        m=self.module(); cells={(1,2,3),(2,2,3)}
        v,f=m.boundary_mesh(cells,depth=0,origin=(0.,0.,0.),side=1.)
        from heart_top import integrate
        expected=[2.,4.,5.,7.,26/3,38/3,74/3,17.5,14.,10.]
        np.testing.assert_allclose(integrate(v,f),expected,atol=1e-12,rtol=0)
        np.testing.assert_allclose(m.box_moments(cells,depth=0,origin=(0.,0.,0.),side=1.),expected,atol=1e-12,rtol=0)
        from make_it_stand.verification.spin_common import native
        facts=native().make_exact_material_region(v,f).facts
        self.assertTrue(facts['completed'] and facts['valid'],facts)

    def test_repaired_material_constraints_are_rechecked_not_copied(self):
        self.assertIsNotNone(importlib.util.find_spec('seal_retention'),'missing final material recheck')
        from seal_retention import top_constraints
        # Aligned asymmetric inertia, COM on +z. Literal geometric moments.
        good=np.array([2.,0.,0.,3.,1.,2.,5.,0.,0.,0.])
        self.assertLess(top_constraints(good,phi=0.,root_side=4.)['maximum_normalized_residual'],1e-12)
        for index in [1,2,7,8,9]:
            bad=good.copy(); bad[index]=.1
            with self.subTest(index=index),self.assertRaisesRegex(ValueError,'CONSTRAINT'):
                top_constraints(bad,phi=0.,root_side=4.)
        bad=good.copy(); bad[1]=np.nan
        with self.assertRaises(ValueError): top_constraints(bad,phi=0.,root_side=4.)


if __name__=='__main__': unittest.main()
