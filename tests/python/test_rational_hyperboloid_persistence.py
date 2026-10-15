"""Actual hyperboloid map/complement/reload gate, not a new optimization run."""
from dataclasses import asdict
from fractions import Fraction
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from make_it_stand.spin import (SpinEvaluatorConfig, SpinItFasterRequest, _native_request,
    evaluate_spin_it_faster_fullrank_quadric)
from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
from make_it_stand.verification.rational_text import fraction_from_rational_text
from make_it_stand.verification.spin_problem import MeshSnapshot
from make_it_stand.verification.spin_bundle import _Writer, _Reader, SCHEMA
from make_it_stand.verification.spin_common import canonical
from test_spin2014_quadric_interface import box


def _new_evidence_directory(root):
    # CTest starts Python in the build directory, not its controller's job
    # directory. Never assume that cwd is exclusive; keep each run durable.
    return Path(tempfile.mkdtemp(prefix='hyperboloid-authority-reload-', dir=root))


class RationalHyperboloidPersistence(unittest.TestCase):
    def test_evidence_directory_is_unique_without_overwriting_prior_payload(self):
        with tempfile.TemporaryDirectory() as root:
            first = _new_evidence_directory(Path(root))
            payload = first / 'problem.json'
            with payload.open('x') as stream:
                stream.write('retained-prior-evidence')
            second = _new_evidence_directory(Path(root))
            self.assertNotEqual(first, second)
            self.assertEqual(first.parent, Path(root))
            self.assertEqual(second.parent, Path(root))
            self.assertEqual(payload.read_text(), 'retained-prior-evidence')
            self.assertEqual(list(second.iterdir()), [])

    def test_actual_production_source_complement_and_saved_reload(self):
        # Same frozen local state as the native rational hyperboloid gate.
        # Its 180 s core test is unchanged. This separate 600 s integration
        # budget also includes exact complement, partition and persistence.
        outer, faces = box(-2., 2.)
        design, df = box(-1., 1.)
        request = SpinItFasterRequest(outer, faces, design, df,
                                      material_density=1., minimum_wall_thickness=.5)
        config = SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1, truncation_radius=8., rational_source_coordinates=True)
        q = (-.2, .07, -.04, .03, 1.2, .9, -1.1, .02, .03, -.04)
        print('[GATE] actual hyperboloid production map begins', flush=True)
        value = evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=config)
        # Save the actual returned map BEFORE independent checks. A downstream
        # failure must retain this expensive producer output for read-only
        # diagnosis; no manual pass flag or optimizer status is attached.
        directory = _new_evidence_directory(Path.cwd())
        print('[EVIDENCE] '+str(directory), flush=True)
        evidence = dict(scope='PRODUCTION_MAP_COMPLEMENT_RELOAD_NOT_OPTIMIZATION',
            request=asdict(request), evaluator=config.identity_dict(), coefficients=q,
            evaluation=asdict(value))
        fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            source = _Writer(fd).write('problem.json', canonical(evidence).encode('ascii'))
        finally:
            os.close(fd)
        self.assertTrue(value.derivative_defined, value.derivative_diagnostic)
        self.assertEqual(np.asarray(value.coefficient_jacobian).shape, (10, 10))
        self.assertTrue(np.isfinite(value.coefficient_jacobian).all())
        authority = value.rational_authority
        for row, lo, hi in zip(authority['source_vertices'],
                authority['source_coordinate_lower'], authority['source_coordinate_upper']):
            for x, lower, upper in zip(row, lo, hi):
                self.assertLessEqual(fraction_from_rational_text(lower), fraction_from_rational_text(x))
                self.assertLessEqual(fraction_from_rational_text(x), fraction_from_rational_text(upper))
                self.assertLessEqual(fraction_from_rational_text(upper)-fraction_from_rational_text(lower),
                                     Fraction(1, 2**384))
        actual = integrate_exact_rational_raw_moments(authority['material_vertices'],
                                                      value.material_faces)
        self.assertEqual(actual, tuple(fraction_from_rational_text(x)
                                      for x in authority['variable_moments']))
        native, _ = _native_request(request)
        print('[GATE] same retained tape complement begins', flush=True)
        exact = native.libigl_exact_boolean_construction_from_rationals(
            authority['design_vertices'], df, authority['material_vertices'],
            value.material_faces, 'difference')
        self.assertTrue(exact['succeeded'])
        retained = MeshSnapshot.from_rationals(authority['material_vertices'], value.material_faces)
        cavity = MeshSnapshot.from_rationals(exact['vertices'], exact['faces'])
        whole = MeshSnapshot.from_rationals(authority['design_vertices'], df)
        partition = native.check_exact_material_partition(whole.exact_region(),
                                                          [retained.exact_region(), cavity.exact_region()])
        self.assertTrue(partition['completed'] and partition['complete_partition'], partition)
        cavity_moments = integrate_exact_rational_raw_moments(cavity.rational_vertices,
                                                              cavity.arrays()[1].tolist())
        whole_moments = integrate_exact_rational_raw_moments(whole.rational_vertices, df)
        self.assertEqual(tuple(a+b for a, b in zip(actual, cavity_moments)), whole_moments)
        # Durable evidence inside this guarded job's unique directory; no
        # optimization status is fabricated, no previous bundle overwritten.
        fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            writer = _Writer(fd)
            records = [writer.mesh(mesh) for mesh in (whole, retained, cavity)]
            writer.write('manifest.json', canonical(dict(source=source, meshes=records,
                                                        partition=partition)).encode('ascii'))
            reader = _Reader(fd)
            self.assertEqual(json.loads(reader.checked(source)), json.loads(canonical(evidence)))
            loaded = [reader.mesh(record, schema=SCHEMA) for record in records]
            self.assertEqual([x.identity for x in loaded],
                             [x.identity for x in (whole, retained, cavity)])
            self.assertEqual(integrate_exact_rational_raw_moments(
                loaded[1].rational_vertices, loaded[1].arrays()[1].tolist()), actual)
        finally:
            os.close(fd)
        print('[GATE] hyperboloid exact authority/complement/reload passed', flush=True)


if __name__ == '__main__':
    unittest.main()
