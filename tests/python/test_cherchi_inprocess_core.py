"""Production in-process Cherchi: actual source, Boolean and derivatives."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import unittest
from fractions import Fraction

import numpy as np
from make_it_stand import spin

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / 'experiments/faster_core_repair_20261012'


class InProcessCoreTests(unittest.TestCase):
    def evaluator(self):
        name = 'make_it_stand.verification.cherchi_inprocess'
        self.assertIsNotNone(importlib.util.find_spec(name),
                             'formal in-process Cherchi evaluator is missing')
        from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
        config = json.loads((CORE / 'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())
        payload = dict(config['request'])
        payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
        return InProcessCherchiEvaluator(spin.SpinItFasterRequest(**payload),
                                       spin.SpinEvaluatorConfig(**config['evaluator']))

    def coefficients(self):
        result = json.loads((CORE / 'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/numerical-result.json').read_text())
        return np.asarray(result['result']['quadric_coefficients'])

    def test_caller_cannot_mutate_bound_physical_design(self):
        original=self.evaluator()
        from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
        physical=np.array(original.request.design_vertices)
        expected=physical.copy()
        evaluator=InProcessCherchiEvaluator(original.request,original.settings,design_vertices_physical=physical)
        physical[0,0]+=1
        np.testing.assert_array_equal(evaluator.design,expected)
        with self.assertRaises(ValueError):evaluator.design[0,0]+=1

    def test_actual_map_preserves_same_call_material_and_complement_source(self):
        evaluator = self.evaluator()
        q = self.coefficients()  # Fixed regression only, never a solver initializer.
        mapped = evaluator(q)
        material = evaluator.material(q)
        self.assertEqual(material['geometry_hash'], mapped['geometry_hash'])
        self.assertEqual(material['quadric_coefficients'], tuple(q))
        self.assertEqual(np.asarray(material['moment_jacobian']).shape, (10, 10))
        self.assertLess(np.max(np.abs(np.asarray(mapped['raw_moments']) -
                                     np.asarray(material['variable_raw_moments']) - evaluator.fixed)), 1e-12)
        self.assertEqual(evaluator.call_counts['boolean'], 1)
        self.assertEqual(evaluator.call_counts['source'], 1)
        self.assertTrue(material['boolean_backend_identity'].startswith('cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:'))
        complement = evaluator.complement(material)
        self.assertEqual(complement['quadric_source_identity'], material['quadric_source_identity'])
        np.testing.assert_array_equal(complement['quadric_source_vertices'], material['quadric_source_vertices'])
        self.assertNotEqual(complement['boolean_operation'], material['boolean_operation'])
        self.assertEqual(evaluator.call_counts['source'], 1)
        self.assertEqual(evaluator.call_counts['boolean'], 2)
        # Returned views cannot change the cached accepted material.
        material['variable_vertices'][0][0] += 1
        self.assertEqual(evaluator.material(q)['geometry_hash'], mapped['geometry_hash'])
        self.assertNotEqual(evaluator.material(q)['variable_vertices'][0][0], material['variable_vertices'][0][0])

    def test_mismatched_coefficients_and_terminal_source_are_rejected(self):
        evaluator = self.evaluator()
        q = self.coefficients()
        evaluator(q)
        with self.assertRaisesRegex(ValueError, 'MATERIAL_NOT_EVALUATED'):
            evaluator.material(q + 0.001)
        material = evaluator.material(q)
        material['quadric_source_vertices'][0][0] += 1
        with self.assertRaisesRegex(ValueError, 'SOURCE_IDENTITY'):
            evaluator.complement(material)

    def test_feature_identity_ignores_output_renumbering_and_continuous_parameters(self):
        self.evaluator()
        from make_it_stand.verification.cherchi_inprocess import feature_branch_identity
        records = [dict(point_type='EXPLICIT3D', parameter_t=0., operands=[dict(
            status='UNIQUE', source_vertices=[dict(input_object=0,input_vertex=i)])]) for i in range(3)]
        a = feature_branch_identity('intersection', [[0,1,2]], [[0,1,2]], 'ELLIPSOID', 0, [], records)
        permuted = [records[2],records[0],records[1]]
        permuted[0]['parameter_t'] = 0.9
        b = feature_branch_identity('intersection', [[1,2,0]], [[0,1,2]], 'ELLIPSOID', 0, [], permuted)
        self.assertEqual(a,b)
        with self.assertRaisesRegex(ValueError, 'CORRESPONDENCE'):
            broken = json.loads(json.dumps(records))
            broken[0]['operands'][0]['status'] = 'AMBIGUOUS'
            feature_branch_identity('intersection', [[0,1,2]], [[0,1,2]], 'ELLIPSOID', 0, [], broken)

    def test_nonidentity_frame_material_returns_to_original_model_coordinates(self):
        original=self.evaluator()
        from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
        axes=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
        evaluator=InProcessCherchiEvaluator(replace(original.request,target_axes=axes.tolist()),original.settings)
        q=self.coefficients()
        evaluator(q)
        material=evaluator.material(q)
        target=(np.asarray(material['variable_vertices_model'])-evaluator.request.target_origin)@axes
        np.testing.assert_allclose(target,material['variable_vertices'],atol=1e-14,rtol=0)

    def test_translated_target_frame_does_not_change_exact_design_partition(self):
        original=self.evaluator()
        from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
        from make_it_stand.verification.spin_problem import MeshSnapshot, partition
        request=replace(original.request,
            outer_vertices=(np.asarray(original.request.outer_vertices)/2.2).tolist(),
            design_vertices=(np.asarray(original.request.design_vertices)/2.2).tolist(),
            target_origin=(np.asarray(original.request.target_origin)/2.2).tolist())
        evaluator=InProcessCherchiEvaluator(request,original.settings)
        q=self.coefficients();evaluator(q)
        material=evaluator.material(q);complement=evaluator.complement(material)
        def mesh(record):
            # Independent exact inverse of the identity-axis target translation.
            rows=[[str(Fraction(x)+Fraction(float(request.target_origin[i])))
                   for i,x in enumerate(row)] for row in record['variable_rational_vertices']]
            return MeshSnapshot.from_rationals(rows,record['variable_faces'])
        domain=MeshSnapshot((request.design_vertices,request.design_faces))
        verdict=partition(domain.exact_region(),[mesh(material).exact_region(),mesh(complement).exact_region()])
        self.assertTrue(verdict['complete_partition'])

    def test_nonbinary_normalization_preserves_physical_input_partition(self):
        original=self.evaluator()
        from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
        from make_it_stand.verification.spin_problem import MeshSnapshot, partition
        request=replace(original.request,
            outer_vertices=(np.asarray(original.request.outer_vertices)/2.2).tolist(),
            design_vertices=(np.asarray(original.request.design_vertices)/2.2).tolist(),
            target_origin=(np.asarray(original.request.target_origin)/2.2).tolist())
        evaluator=InProcessCherchiEvaluator(request,original.settings,
            normalization=dict(origin=[0.,0.,0.],L=2.2),
            design_vertices_physical=original.request.design_vertices)
        q=self.coefficients();evaluator(q)
        retained=evaluator.material(q);carved=evaluator.complement(retained)
        def mesh(record):
            def text(value):return f'{value.numerator}/{value.denominator}'
            rows=[[text((Fraction(x)+Fraction(float(request.target_origin[i])))*Fraction(2.2))
                   for i,x in enumerate(row)] for row in record['variable_rational_vertices']]
            return MeshSnapshot.from_rationals(rows,record['variable_faces'])
        domain=MeshSnapshot((original.request.design_vertices,request.design_faces))
        self.assertTrue(partition(domain.exact_region(),[mesh(retained).exact_region(),mesh(carved).exact_region()])['complete_partition'])


if __name__ == '__main__':
    unittest.main()
