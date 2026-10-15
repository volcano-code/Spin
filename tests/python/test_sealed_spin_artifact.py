"""Real exact-native analytic fixtures; these do not admit original paper assets."""
import importlib
from fractions import Fraction as Q
import unittest

import numpy as np

from test_exact_material_region import box, combine, l_prism


ORDER = ['M', 'LX', 'LY', 'LZ', 'QXX', 'QYY', 'QZZ', 'QYZ', 'QXZ', 'QXY']
TARGET = {'mode': 'SPINNING_TOP', 'axis': [0., 0., 1.], 'origin': [0., 0., 0.],
          'support_points': []}


class SpinSealTests(unittest.TestCase):
    def setUp(self):
        name = 'make_it_stand.verification.spin_artifact'
        self.assertIsNotNone(importlib.util.find_spec(name), 'integrated spin seal is missing')
        self.a = importlib.import_module(name)
        self.n = importlib.import_module('_make_it_stand')

    def problem(self, outer=None, fixed=(), design=None, wall=.25, **kw):
        outer = box() if outer is None else outer
        return self.a.SpinArtifactProblem(original=outer, fixed_regions=fixed,
            design_regions=[('d', outer, 1.)] if design is None else design,
            length_unit='m', mass_unit='kg', target=TARGET, minimum_outer_wall=wall, **kw)

    def claim(self, p, regions):
        # Analytic rectangular-solid products, not either surface integrator.
        # L = [0,3]x[0,1]x[0,1] union [0,1]x[1,3]x[0,1].
        # The remaining fixtures consist of ordered eight-vertex boxes.
        frame = p.normalization
        total = np.zeros(10)
        published = np.zeros(10)
        have_published = hasattr(self.n, 'differentiate_surface_raw_moments')
        for mesh, density in regions:
            if mesh is not None:
                v, f = mesh
                if len(v) == 12:
                    boxes = [((0,0,0), (3,1,1), 1), ((0,1,0), (1,3,1), 1)]
                else:
                    self.assertEqual(len(v) % 8, 0)
                    boxes = []
                    for i in range(len(v)//8):
                        vertices = v[i*8:(i+1)*8]
                        sign = 1 if list(f[i*12]-8*i) == [0,2,1] else -1
                        boxes.append((vertices.min(axis=0), vertices.max(axis=0), sign))
                exact = [Q(0) for _ in range(10)]
                for lo, hi, sign in boxes:
                    low = [(Q(float(x))-Q(float(o)))/Q(frame['L']) for x,o in zip(lo, frame['origin'])]
                    high = [(Q(float(x))-Q(float(o)))/Q(frame['L']) for x,o in zip(hi, frame['origin'])]
                    volume = Q(sign)*Q(float(density))/Q(frame['R'])
                    for a,b in zip(low,high): volume *= b-a
                    mean = [(a+b)/2 for a,b in zip(low,high)]
                    diag = [(a*a+a*b+b*b)/3 for a,b in zip(low,high)]
                    terms = [volume, *[volume*x for x in mean], *[volume*x for x in diag],
                             volume*mean[1]*mean[2], volume*mean[0]*mean[2], volume*mean[0]*mean[1]]
                    exact = [a+b for a,b in zip(exact,terms)]
                total += np.array([float(x) for x in exact])
                if have_published:
                    published += self.n.differentiate_surface_raw_moments(
                        (v-np.array(frame['origin']))/frame['L'], f, density/frame['R'])['moments']
        if have_published:
            np.testing.assert_allclose(published, total, atol=1e-10, rtol=0)
            total = published  # Separate actual native claim, never checker values.
        return {'order': ORDER, 'frame': 'INPUT_NORMALIZED_V1', 'moments': total.tolist()}

    def seal(self, p, outer, retained, carved, material, carved_output, claim=None):
        return self.a.SealedSpinArtifact.create(p, outer=outer,
            design_results=[('d', retained, carved)], final_material=material,
            newly_carved=carved_output,
            solver_claim=self.claim(p, [(retained, 1.)]) if claim is None else claim)

    def hollow(self, wall=.25):
        outer, cavity = box(), box((.25,)*3, (.75,)*3)
        retained = combine(outer, (cavity[0], cavity[1][:, ::-1]))
        p = self.problem(outer, wall=wall)
        return p, outer, retained, cavity

    def test_hollow_box_exact_wall_and_independent_literal_moments(self):
        p, o, r, c = self.hollow()
        s = self.seal(p, o, r, c, r, c)
        f = s.facts
        self.assertEqual(f['geometry_status'], 'MESH_MATERIAL_VALIDATED')
        self.assertEqual(f['validation_scope'], 'HOLLOW_ONLY_MESH_MATERIAL_V1')
        self.assertEqual(f['clearance']['squared_distance_exact'], '1/16')
        self.assertEqual(f['wall_thickness_audit'], {
            'outer_to_newly_carved': 'VERIFIED',
            'same_component_feature': 'UNKNOWN',
            'inter_cavity_web': 'NOT_APPLICABLE',
            'overall': 'UNKNOWN'})
        self.assertEqual(f['manufacturing_admission'], 'BLOCKED_UNKNOWN_WALL_THICKNESS')
        np.testing.assert_allclose(f['normalized']['moments'],
            [7/8, 0, 0, 0, 31/384, 31/384, 31/384, 0, 0, 0], atol=1e-10, rtol=0)
        self.assertAlmostEqual(f['physical']['mass'], 7/8)
        np.testing.assert_allclose(f['physical']['center_of_mass'], [.5]*3, atol=1e-10)
        self.assertEqual((f['material_components'], f['bounded_void_components'],
                          f['newly_carved_components']), (1, 1, 1))
        self.assertEqual(f['solver_target_feasibility'], 'UNKNOWN')
        self.assertEqual(f['build_source_attestation'], 'NOT_VERIFIED')
        self.assertEqual(f['physical_validation'], 'NOT_IN_SCOPE')
        self.assertEqual(f['original_evidence_admission'], 'BLOCKED_UNKNOWN')
        self.assertEqual(f['inter_cavity_web_thickness'], 'NOT_APPLICABLE')
        self.assertEqual(f['same_component_feature_thickness'], 'UNKNOWN')
        self.assertEqual(s.revalidate().scientific_digest, s.scientific_digest)

    def test_wall_above_actual_distance_rejected(self):
        p, o, r, c = self.hollow(wall=np.nextafter(.25, 1.))
        with self.assertRaises(ValueError): self.seal(p, o, r, c, r, c)

    def test_no_new_cavity_is_explicit_not_applicable(self):
        o = box(); p = self.problem(o)
        s = self.seal(p, o, o, None, o, None)
        self.assertEqual(s.facts['clearance']['state'], 'NOT_APPLICABLE')
        self.assertIsNone(s.facts['clearance']['squared_distance_exact'])
        self.assertEqual(s.facts['newly_carved_components'], 0)
        self.assertEqual(s.facts['wall_thickness_audit']['outer_to_newly_carved'], 'NOT_APPLICABLE')
        self.assertEqual(s.facts['wall_thickness_audit']['inter_cavity_web'], 'NOT_APPLICABLE')
        self.assertEqual(s.facts['wall_thickness_audit']['overall'], 'UNKNOWN')
        self.assertEqual(s.facts['manufacturing_admission'], 'BLOCKED_UNKNOWN_WALL_THICKNESS')

    def test_nonconvex_l_prism_hand_integrated_properties(self):
        o = l_prism(); p = self.problem(o)
        s = self.seal(p, o, o, None, o, None)
        f = s.facts['physical']
        np.testing.assert_allclose(f['moments'],
            [5, 11/2, 11/2, 5/2, 29/3, 29/3, 5/3, 11/4, 11/4, 17/4], atol=1e-10, rtol=0)
        np.testing.assert_allclose(f['center_of_mass'], [11/10, 11/10, 1/2], atol=1e-10)
        np.testing.assert_allclose(f['central_inertia'],
            [[121/30, 9/5, 0], [9/5, 121/30, 0], [0, 0, 217/30]], atol=1e-10, rtol=0)
        np.testing.assert_allclose(s.facts['normalized']['moments'],
            [5/27, -2/81, -2/81, 0, 53/2916, 53/2916, 5/2916, 0, 0, -1/243], atol=1e-10, rtol=0)

    def test_nested_island_multiple_void_counts_are_actual(self):
        o = box((0,)*3, (8,)*3)
        a, b, c = box((1,)*3, (7,)*3), box((2,)*3, (6,)*3), box((3,)*3, (5,)*3)
        r = combine(o, (a[0], a[1][:, ::-1]), b, (c[0], c[1][:, ::-1]))
        carved = combine(a, (b[0], b[1][:, ::-1]), c)
        p = self.problem(o, wall=1.)
        s = self.seal(p, o, r, carved, r, carved)
        self.assertEqual((s.facts['material_components'], s.facts['bounded_void_components'],
                          s.facts['newly_carved_components']), (2, 2, 2))
        self.assertEqual(s.facts['wall_thickness_audit']['inter_cavity_web'], 'UNKNOWN')

    def test_two_touching_densities_and_fixed_material(self):
        o = box(); left = box(hi=(.5, 1, 1)); right = box(lo=(.5, 0, 0))
        p = self.problem(o, fixed=[('f', left, 2.)], design=[('d', right, 4.)])
        claim = self.claim(p, [(left, 2.), (right, 4.)])
        s = self.seal(p, o, right, None, o, None, claim)
        self.assertAlmostEqual(s.facts['physical']['mass'], 3.)
        np.testing.assert_allclose(s.facts['physical']['center_of_mass'], [7/12, .5, .5], atol=1e-10)
        self.assertEqual([x['density'] for x in s.facts['regions']], [2., 4.])
        self.assertEqual(p.normalization['R'], 4.)
        with self.assertRaises(ValueError): self.seal(p, o, right, None, right, None, claim)
        wrong = self.claim(p, [(left, 4.), (right, 16.)])
        with self.assertRaises(ValueError): self.seal(p, o, right, None, o, None, wrong)

    def test_problem_rejects_missing_leaking_overlapping_partitions(self):
        o = box()
        for parts in [[('d', box(hi=(.5, 1, 1)), 1.)],
                      [('d', box(hi=(2, 1, 1)), 1.)],
                      [('d', o, 1.), ('e', o, 1.)], []]:
            with self.subTest(parts=len(parts)), self.assertRaises(ValueError):
                self.problem(o, design=parts)

    def test_final_partition_missing_overlap_outside_and_combined_mismatch(self):
        p, o, r, c = self.hollow()
        for retained, carved, material, combined in [
            (r, None, r, None), (o, c, o, c),
            (o, box((2,)*3, (3,)*3), o, box((2,)*3, (3,)*3)),
            (r, c, o, c), (r, c, r, None), (None, o, None, o)]:
            with self.subTest(), self.assertRaises(ValueError):
                self.seal(p, o, retained, carved, material, combined,
                          self.claim(p, [(o, 1.)]))

    def test_original_signed_zero_and_order_bytes_must_match(self):
        p, o, r, c = self.hollow()
        changed = (o[0].copy(), o[1].copy()); changed[0][0, 0] = -0.
        for outer in [changed, (o[0], o[1][::-1])]:
            with self.assertRaises(ValueError): self.seal(p, outer, r, c, r, c)

    def test_original_preexisting_void_and_multiple_components(self):
        o = combine(box(), box((.25,)*3, (.75,)*3, True), box((2,0,0), (3,1,1)))
        p = self.problem(o)
        s = self.seal(p, o, o, None, o, None)
        self.assertEqual((s.facts['material_components'], s.facts['bounded_void_components']), (2, 1))

    def test_invalid_open_crossing_winding_and_parser_loss_rejected(self):
        o = box(); integer = o[0].astype(np.int64); integer[1, 0] = 2**53+1
        cases = [(o[0], o[1][:-1]), combine(o, box((.5,0,0),(1.5,1,1))),
                 (o[0], o[1][:, ::-1]), (o[0], o[1].astype(float)),
                 (o[0].astype(bool), o[1]), (integer, o[1]),
                 (np.full_like(o[0], np.nan), o[1]),
                 (np.empty((0,3)), np.empty((0,3), dtype=int))]
        if np.dtype(np.longdouble).itemsize > 8: cases.append((o[0].astype(np.longdouble), o[1]))
        for mesh in cases:
            with self.subTest(), self.assertRaises((ValueError, TypeError)):
                self.problem(mesh)

    def test_all_ten_wrong_claims_order_frame_nonfinite_rejected(self):
        p, o, r, c = self.hollow()
        good = self.claim(p, [(r, 1.)])
        for i in range(10):
            wrong = dict(good, moments=good['moments'].copy()); wrong['moments'][i] += 2e-6
            with self.subTest(i=i), self.assertRaises(ValueError): self.seal(p, o, r, c, r, c, wrong)
        for wrong in [dict(good, frame='WORLD'), dict(good, order=ORDER[::-1]),
                      dict(good, moments=[float('nan')]*10), dict(good, moments=[float('inf')]*10),
                      dict(good, moments=[-1.]+[0.]*9), dict(good, passed=True)]:
            with self.subTest(), self.assertRaises(ValueError): self.seal(p, o, r, c, r, c, wrong)

    def test_missing_duplicate_unknown_final_ids_and_density_override(self):
        p, o, r, c = self.hollow()
        for results in [[], [('d',r,c), ('d',r,c)], [('x',r,c)], [('d',r,c,2.)]]:
            with self.subTest(), self.assertRaises(ValueError):
                self.a.SealedSpinArtifact.create(p, outer=o, design_results=results,
                    final_material=r, newly_carved=c, solver_claim=self.claim(p, [(r,1.)]))

    def test_problem_metadata_and_unsupported_profiles_fail(self):
        for wall in [0., -1., float('nan'), float('inf'), True]:
            with self.subTest(wall=wall), self.assertRaises((ValueError, TypeError)): self.problem(wall=wall)
        for density in [0., -1., float('inf'), True]:
            with self.subTest(), self.assertRaises((ValueError, TypeError)):
                self.problem(design=[('d', box(), density)])
        for profile in ['CAD', 'HYDROSTATIC', 'DEFORMATION']:
            with self.subTest(), self.assertRaises(ValueError): self.problem(validation_scope=profile)

    def test_input_and_facts_copies_cannot_mutate_seal(self):
        p, o, r, c = self.hollow(); digest = p.input_digest
        s = self.seal(p, o, r, c, r, c); before = s.scientific_digest
        o[0][:] = 42.; p.original.arrays()[0][:] = 99.
        p.normalization['origin'][0] = 88.
        s.facts['physical']['mass'] = -1.
        self.assertEqual(p.input_digest, digest)
        self.assertEqual(s.revalidate().scientific_digest, before)
        with self.assertRaises(AttributeError): p.minimum_outer_wall = 0.
        with self.assertRaises(AttributeError): s.scientific_digest = 'fake'

    def test_target_modes_support_data_units_and_normalization(self):
        o = box()
        kwargs = dict(original=o, fixed_regions=[], design_regions=[('d',o,1.)],
                      length_unit='mm', mass_unit='g', minimum_outer_wall=.25)
        for target in [dict(TARGET, mode='unknown'), dict(TARGET, axis=[0.,0.,0.]),
                       dict(TARGET, origin=[float('inf'),0.,0.]),
                       dict(TARGET, mode='STATIC_BALANCE'), dict(TARGET, support_points=[[0.,0.]]),
                       dict(TARGET, passed=True)]:
            with self.subTest(), self.assertRaises(ValueError): self.a.SpinArtifactProblem(**kwargs, target=target)
        target = dict(TARGET, mode='STATIC_BALANCE', support_points=[[0.,0.,0.]], axis=[0.,0.,2.])
        p = self.a.SpinArtifactProblem(**kwargs, target=target)
        self.assertEqual(p.target['axis'], [0.,0.,2.])
        self.assertEqual(p.identity['fixed_regions'], [])
        for key in ['length_unit', 'mass_unit']:
            bad = dict(kwargs); bad[key] = ''
            with self.assertRaises(ValueError): self.a.SpinArtifactProblem(**bad, target=TARGET)

    def test_invalid_final_mesh_and_positive_mass_requirement(self):
        p, o, r, c = self.hollow()
        with self.assertRaises(ValueError): self.seal(p, o, (r[0],r[1][:-1]), c, r, c)
        with self.assertRaises(ValueError): self.seal(p, o, r, c, (r[0],r[1][:-1]), c)
        with self.assertRaises(ValueError): self.seal(p, o, None, o, None, o, self.claim(p,[(o,1.)]))

    def test_original_scale_large_origin_and_normalized_density_underflow(self):
        o = box((2**40,)*3, (2**40+1,)*3)
        p = self.problem(o)
        s = self.seal(p,o,o,None,o,None)
        np.testing.assert_allclose(s.facts['normalized']['moments'],
                                   [1,0,0,0,1/12,1/12,1/12,0,0,0], atol=1e-10, rtol=0)
        np.testing.assert_allclose(s.facts['physical']['central_inertia'], np.eye(3)/6, atol=1e-10, rtol=0)
        for scale in [1e-200, 1e100]:
            tiny = (box()[0]*scale, box()[1])
            with self.assertRaises(ValueError): self.problem(tiny)
        left, right = box(hi=(.5,1,1)), box(lo=(.5,0,0))
        with self.assertRaises(ValueError): self.problem(design=[('d',left,5e-324), ('e',right,1e308)])

    def test_public_problem_artifact_bundle_api(self):
        public = importlib.import_module('make_it_stand.verification')
        for name in ('SpinArtifactProblem', 'SealedSpinArtifact', 'save_spin_artifact', 'load_spin_artifact'):
            self.assertTrue(hasattr(public, name), 'missing public API: ' + name)

    def test_two_design_ids_and_explicit_empty_retained(self):
        o = box(); c = box((.25,)*3,(.75,)*3)
        r = combine(o, (c[0], c[1][:,::-1]))
        p = self.problem(o, design=[('shell',r,2.), ('core',c,7.)])
        s = self.a.SealedSpinArtifact.create(p, outer=o,
            design_results=[('shell',r,None), ('core',None,c)], final_material=r, newly_carved=c,
            solver_claim=self.claim(p,[(r,2.)]))
        self.assertAlmostEqual(s.facts['physical']['mass'], 7/4)
        self.assertTrue(s.facts['regions'][1]['explicit_empty'])
        self.assertEqual(p.normalization['R'], 7.)
        for bad in [[('shell',r,None)], [('shell',r,None),('core',None,c),('extra',None,None)]]:
            with self.assertRaises(ValueError):
                self.a.SealedSpinArtifact.create(p, outer=o, design_results=bad,
                    final_material=r, newly_carved=c, solver_claim=self.claim(p,[(r,2.)]))


class SpinSequenceIngestionTests(unittest.TestCase):
    """Inspect supplied scalar types before NumPy can erase their provenance."""
    def setUp(self):
        from make_it_stand.verification.spin_problem import MeshSnapshot
        self.snapshot = MeshSnapshot

    def assert_snapshot_bytes(self, mesh, expected):
        try:
            saved = self.snapshot(mesh).arrays()
        except (TypeError, ValueError) as exc:
            self.fail('representable mixed mesh was rejected: ' + str(exc))
        self.assertEqual(saved[0].tobytes(), np.asarray(expected[0], dtype='<f8').tobytes())
        self.assertEqual(saved[1].tobytes(), np.asarray(expected[1], dtype='<i4').tobytes())

    def test_mixed_list_tuple_python_numpy_bools_rejected(self):
        for boolean in (False, True, np.bool_(False), np.bool_(True)):
            for sequence in (list, tuple):
                for field in (0, 1):
                    v, f = (a.tolist() for a in box())
                    if field == 0:
                        v[1 if boolean else 0][0] = boolean
                    else:
                        f[0][2 if boolean else 0] = boolean
                    mesh = sequence((sequence(sequence(row) for row in v),
                                     sequence(sequence(row) for row in f)))
                    with self.subTest(type=type(boolean).__name__, sequence=sequence.__name__, field=field), \
                         self.assertRaises((ValueError, TypeError)):
                        self.snapshot(mesh)

    def test_lossy_mixed_python_numpy_integers_rejected(self):
        for lower in (9007199254740993, np.int64(9007199254740993),
                      np.uint64(9007199254740993), -9007199254740993):
            for sequence in (list, tuple):
                v, f = box()
                upper = float(int(lower) + 4)
                vertices = sequence(sequence((lower if x == 0 else upper, float(y), float(z)))
                                    for x, y, z in v)
                with self.subTest(lower=lower, sequence=sequence.__name__), \
                     self.assertRaises((ValueError, TypeError)):
                    self.snapshot((vertices, f.tolist()))

    def test_representable_mixed_rows_and_integer_index_types_preserved(self):
        v, f = box(); v[0, 0] = -0.
        vertices = v.tolist()
        vertices[0] = (np.float32(-0.), np.int64(0), 0.)
        vertices[1] = [1, np.float16(0.), np.uint32(0)]
        vertices[2] = np.array(vertices[2], dtype=np.float32)
        faces = f.tolist()
        faces[0] = (np.uint64(0), np.int64(2), 1)
        faces[1] = np.array(faces[1], dtype=np.int16)
        for sequence in (list, tuple):
            self.assert_snapshot_bytes(sequence((sequence(vertices), sequence(faces))), (v, f))
        self.assertIsInstance(vertices[0][0], np.float32)
        self.assertTrue(np.signbit(vertices[0][0]))
        self.assertIsInstance(faces[0][0], np.uint64)

    def test_large_exactly_representable_integer_float_mixtures_preserved(self):
        v, f = box((9007199254740994, 0, 0), (9007199254740998, 1, 1))
        vertices = [(np.uint64(9007199254740994) if x == 9007199254740994 else float(x),
                     np.float32(y), int(z)) for x, y, z in v]
        self.assert_snapshot_bytes((vertices, f.tolist()), (v, f))

    def test_unsupported_sequence_scalars_precision_and_index_bounds_rejected(self):
        from decimal import Decimal
        bad_vertices = [Decimal('0'), Q(0), '0', 0j, None, float('nan'), float('inf')]
        if np.dtype(np.longdouble).itemsize > 8:
            bad_vertices.append(np.longdouble(0))
        for value in bad_vertices:
            v, f = (a.tolist() for a in box()); v[0][0] = value
            with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                self.snapshot((v, f))
        for value in (0., np.float32(0), -1, 2**31, np.uint64(2**64-1), '0', None):
            v, f = (a.tolist() for a in box()); f[0][0] = value
            with self.subTest(index=value), self.assertRaises((ValueError, TypeError)):
                self.snapshot((v, f))

    def test_direct_arrays_and_noncontiguous_views_keep_native_validation(self):
        v, f = box(); v[0, 0] = -0.
        for dtype in (np.float16, np.float32, np.float64):
            coordinates = np.asfortranarray(v.astype(dtype))
            self.assert_snapshot_bytes((coordinates, f[:, ::-1][:, ::-1]), (v, f))
        for dtype in (np.int16, np.int64, np.uint64):
            coordinates = box()[0].astype(dtype)
            self.assert_snapshot_bytes((coordinates, f.astype(dtype)), box())
        for bad in ((v.astype(bool), f), (v.astype(object), f), (v, f.astype(float)),
                    (v, f.astype(bool)), (v, f.astype(object))):
            with self.subTest(), self.assertRaises((ValueError, TypeError)):
                self.snapshot(bad)

    def test_problem_and_final_output_reject_erased_boolean_type(self):
        fixture = SpinSealTests(); fixture.setUp()
        v, f = box(); coordinates = v.tolist(); coordinates[0][0] = False
        with self.assertRaises((ValueError, TypeError)):
            fixture.problem((coordinates, f.tolist()))
        p = fixture.problem((v, f))
        with self.assertRaises((ValueError, TypeError)):
            fixture.seal(p, (coordinates, f.tolist()), (v, f), None, (v, f), None)


if __name__ == '__main__': unittest.main()
