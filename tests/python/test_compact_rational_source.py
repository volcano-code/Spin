"""Actual opt-in material mapping; no optimizer or FD completion claim."""
from dataclasses import replace
from fractions import Fraction as F
import json
import unittest

from make_it_stand.spin import (SpinEvaluatorConfig, SpinItFasterRequest,
    evaluate_spin_it_faster_fullrank_quadric)
from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
from test_spin2014_quadric_interface import box


class CompactRationalSource(unittest.TestCase):
    def test_policy_is_explicit_and_rejects_invalid_dependencies(self):
        # Catch unbound new policy, unsafe defaults and compact-with-float use.
        old=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
                               rational_source_coordinates=True)
        self.assertFalse(old.compact_rational_source_coordinates)
        self.assertNotIn('compact_rational_source_coordinates',old.identity_dict())
        new=replace(old,compact_rational_source_coordinates=True)
        self.assertTrue(new.identity_dict()['compact_rational_source_coordinates'])
        with self.assertRaises(ValueError):
            SpinEvaluatorConfig(compact_rational_source_coordinates=True)
        with self.assertRaises(ValueError):
            replace(old,compact_rational_source_coordinates=1)

    def test_real_mapping_preserves_bounds_and_binds_alternate_authority(self):
        # Catch binding omission, floating reconstitution or stale source hash.
        outer,faces=box(-2.,2.)
        design,df=box(-1.,1.)
        request=SpinItFasterRequest(outer,faces,design,df,material_density=2.7)
        old=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=0,rational_source_coordinates=True)
        raw=(-.25,0.,0.,0.,1.,1.,1.,0.,0.,0.)
        baseline=evaluate_spin_it_faster_fullrank_quadric(request,raw,evaluator_settings=old)
        compact=evaluate_spin_it_faster_fullrank_quadric(request,raw,
            evaluator_settings=replace(old,compact_rational_source_coordinates=True))
        a,b=baseline.rational_authority,compact.rational_authority
        self.assertNotEqual(a['source_construction_identity'],b['source_construction_identity'])
        self.assertNotEqual(baseline.geometry_hash,compact.geometry_hash)
        self.assertEqual(a['source_faces'],b['source_faces'])
        for key in ('source_coordinate_lower','source_coordinate_upper'):
            self.assertEqual(a[key],b[key])
        changed=False
        for row,lo,hi,oldrow in zip(b['source_vertices'],b['source_coordinate_lower'],
                                   b['source_coordinate_upper'],a['source_vertices']):
            for x,l,u,o in zip(row,lo,hi,oldrow):
                x,l,u,o=map(F,(x,l,u,o))
                self.assertLessEqual(l,x);self.assertLessEqual(x,u)
                self.assertEqual((2**1024)%x.denominator,0)
                changed|=x!=o
        self.assertTrue(changed)
        exact=integrate_exact_rational_raw_moments(b['material_vertices'],compact.material_faces)
        self.assertEqual(tuple(F(request.material_density)*x for x in exact),
                         tuple(map(F,b['variable_moments'])))
        self.assertTrue(compact.derivative_defined)
        # Actual stale-policy admission and SAME-retained complement/save.
        from types import SimpleNamespace
        import os
        import tempfile
        from make_it_stand.verification.faster_material import (
            _exact_spin2014_saved_retained,_canonical_exact_boolean_boundary)
        from make_it_stand.verification.spin_bundle import _Writer,_Reader,SCHEMA
        endpoint=SimpleNamespace(quadric_coefficients=compact.normalized_coefficients,
            final_geometry_hash=compact.geometry_hash,material_vertices=compact.material_vertices,
            material_faces=compact.material_faces,rational_authority=b)
        candidate=SimpleNamespace(measured_request=request,solve=endpoint)
        new=replace(old,compact_rational_source_coordinates=True)
        retained,_=_exact_spin2014_saved_retained(candidate,new)
        with self.assertRaises(ValueError):
            _exact_spin2014_saved_retained(candidate,old)
        from make_it_stand.verification.faster_material import native
        from make_it_stand.verification.rational_text import fraction_to_rational_text
        domain=_canonical_exact_boolean_boundary(tuple(tuple(
            fraction_to_rational_text(F(x)) for x in row) for row in request.design_vertices),
            request.design_faces)
        complement=native().libigl_exact_boolean_construction_from_rationals(
            domain.rational_vertices,domain.arrays()[1],
            retained.rational_vertices,retained.arrays()[1],'difference')
        self.assertTrue(complement['succeeded'])
        cavity=_canonical_exact_boolean_boundary(complement['vertices'],complement['faces'])
        with tempfile.TemporaryDirectory() as folder:
            fd=os.open(folder,os.O_RDONLY|os.O_DIRECTORY)
            try:
                writer=_Writer(fd)
                records=[writer.mesh(x) for x in (retained,cavity)]
                reader=_Reader(fd)
                restored=[reader.mesh(x,schema=SCHEMA) for x in records]
                self.assertEqual([x.identity for x in restored],[retained.identity,cavity.identity])
            finally:
                os.close(fd)


if __name__=='__main__':
    unittest.main()
