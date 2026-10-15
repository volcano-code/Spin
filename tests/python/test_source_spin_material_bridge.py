"""Actual source solve, immutable run binding, and independent material seal."""
import importlib
import json
import unittest
from dataclasses import asdict, replace, FrozenInstanceError
import tempfile
from pathlib import Path
import numpy as np
from test_exact_material_region import combine
from make_it_stand.verification.spin_problem import SpinArtifactProblem


def box(lo, hi):
    # Same indexed input as source_surface_test_utils.h, required by native hashes.
    v = np.array([[hi[a] if i & (1 << a) else lo[a] for a in range(3)] for i in range(8)], dtype=float)
    f = np.array([[0,2,1],[1,2,3],[4,5,6],[5,7,6],[0,1,4],[1,5,4],
                  [2,6,3],[3,6,7],[0,4,2],[2,4,6],[1,3,5],[3,7,5]], dtype=np.int32)
    return v, f


def problem(density=1., wall=.5, outer=None, design=None, mode='STATIC_BALANCE', **target):
    o = box((-1,-1,0), (1,1,2)) if outer is None else outer
    d = box((-.5,-.5,1), (.5,.5,1.5)) if design is None else design
    return SpinArtifactProblem(original=o,
        fixed_regions=[('shell', combine(o, (d[0], d[1][:, ::-1])), density)],
        design_regions=[('design', d, density)], length_unit='m', mass_unit='kg',
        target={'mode': mode, 'axis': [0,0,1], 'origin': [0,0,0],
                'support_points': [[0,0,0]] if mode == 'STATIC_BALANCE' else [], **target},
        minimum_outer_wall=wall)


def settings(m):
    return dict(source_settings=m.SourceGeometrySettings(root_origin=(-1,-1,0), root_side=2., maximum_depth=4),
        solver_settings=m.SourceSolverSettings(maximum_optimization_level=4, laplacian_weight=1.,
            slqp=m.AuditedSlqpSettings(maximum_iterations=300, equality_tolerance=1e-10,
                                      projected_gradient_tolerance=1e-8)),
        boundary_settings=m.SourceBoundarySettings())


class SourceBridgeAvailabilityTests(unittest.TestCase):
    def test_actual_source_material_entrypoint_exists(self):
        n = importlib.import_module('_make_it_stand')
        self.assertTrue(hasattr(n, 'solve_source_spin2014_material'),
                        'missing actual native source material solve entrypoint')
        self.assertIsNotNone(importlib.util.find_spec(
                        'make_it_stand.verification.source_spin_material'), 'missing independent sealing bridge')

    def test_readonly_source_domain_inspection_exposes_sparse_laplacian(self):
        """Independent auditors need the actual source-derived leaf graph.

        This is deliberately an inspection-only endpoint: it cannot accept a
        caller-supplied moment matrix or mutate the adaptive partition.
        """
        native = importlib.import_module('_make_it_stand')
        self.assertTrue(hasattr(native, 'inspect_source_spin2014_domain'),
                        'missing read-only source-domain inspection entrypoint')
        outer = box((-1., -1., 0.), (1., 1., 2.))
        design = box((-.5, -.5, 1.), (.5, .5, 1.5))
        observed = native.inspect_source_spin2014_domain(
            *outer, *design,
            source_settings={'root_origin': [-1., -1., 0.], 'root_side': 2.,
                             'maximum_depth': 4, 'maximum_nodes': 5_000_000},
            initial_level=2, maximum_optimization_level=4,
            maximum_leaf_count=1_000_000,
        )
        self.assertEqual(observed['schema'], 'source-spin2014-domain-inspection-v1')
        self.assertGreater(observed['source']['stats']['boundary_nodes'], 0)
        cells = observed['source_cells']
        self.assertEqual(cells[0]['index'], 0)
        self.assertEqual(cells[0]['level'], 0)
        boundary = [cell for cell in cells if cell['kind'] == 'BOUNDARY']
        self.assertEqual(len(boundary), observed['source']['stats']['boundary_nodes'])
        self.assertTrue(all(cell['level'] == 4 for cell in boundary))
        for cell in cells[1:]:
            self.assertLess(cell['parent'], cell['index'])
        for cell in cells:
            if cell['kind'] == 'INTERNAL':
                self.assertEqual(len(cell['children']), 8)
        leaves = observed['leaves']
        self.assertGreater(len(leaves), 1)
        self.assertEqual(observed['leaf_count'], len(leaves))
        self.assertEqual(len(observed['partition_hash']), 64)
        self.assertEqual(len(observed['source']['partition_hash']), 64)
        edge_set = {tuple(edge) for edge in observed['edges']}
        self.assertTrue(edge_set)
        matrix = np.zeros((len(leaves), len(leaves)))
        for row, column, value in observed['laplacian_triplets']:
            matrix[row, column] += value
        np.testing.assert_allclose(matrix, matrix.T, atol=0, rtol=0)
        np.testing.assert_allclose(matrix.sum(axis=1), 0., atol=1e-14, rtol=0)
        for first, second in edge_set:
            self.assertEqual(matrix[first, second], -1.)
            self.assertEqual(matrix[second, first], -1.)
        beta = np.linspace(.1, .9, len(leaves))
        independent_energy = sum((beta[first] - beta[second]) ** 2
                                 for first, second in edge_set)
        self.assertAlmostEqual(float(beta @ matrix @ beta), independent_energy,
                               places=12)

    def test_readonly_source_domain_summary_avoids_expanding_cell_payloads(self):
        """Large source trees can be audited without materializing Python records."""
        native = importlib.import_module('_make_it_stand')
        outer = box((-1., -1., 0.), (1., 1., 2.))
        design = box((-.5, -.5, 1.), (.5, .5, 1.5))
        observed = native.inspect_source_spin2014_domain(
            *outer, *design,
            source_settings={'root_origin': [-1., -1., 0.], 'root_side': 2.,
                             'maximum_depth': 4, 'maximum_nodes': 5_000_000},
            initial_level=2, maximum_optimization_level=4,
            maximum_leaf_count=1_000_000,
            include_details=False,
        )
        self.assertEqual(observed['inspection_detail'], 'SUMMARY')
        self.assertGreater(observed['source_cell_count'], 0)
        self.assertGreater(observed['leaf_count'], 0)
        self.assertTrue(all(type(level) is str
                            for level in observed['leaf_level_counts']))
        self.assertEqual(sum(observed['leaf_level_counts'].values()),
                         observed['leaf_count'])
        self.assertGreater(observed['edge_count'], 0)
        self.assertEqual(len(observed['partition_hash']), 64)
        self.assertEqual(len(observed['source']['source_hash']), 64)
        for expanded in ('source_cells', 'leaves', 'edges', 'laplacian_triplets'):
            self.assertNotIn(expanded, observed)


class SourceMaterialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = importlib.import_module('make_it_stand.verification.source_spin_material')
        cls.n = importlib.import_module('_make_it_stand')
        cls.p = problem()
        cls.cfg = settings(cls.m)
        cls.solved = cls.m.solve_and_seal_source_spin2014(cls.p, **cls.cfg)

    def test_literal_fixture_and_full_trace(self):
        r = self.solved
        self.assertEqual(r.status, 'BINARY_MOMENT_VALIDATED')
        self.assertIsNone(r.seal_error)
        self.assertIsNotNone(r.artifact)
        s = r.identity['native_run']['solve']
        # v2 binds the explicit termination policy and post-transition facts.
        # Scientific trajectory, cuboid count and analytic moments below retain
        # their v1 expectations; this is an identity-format migration only.
        self.assertEqual(s['input_hash'], 'a7c9239fa8717f83d2cef88c838d7741ecc8e658cbee19ace3387520a89f23e0')
        self.assertEqual(s['trace_hash'], 'f4d0604b0c8d76c29e7d24ac850ee9b558e6ec105ebbb0a1193a8d53a119cc87')
        self.assertEqual(s['output_hash'], 'd8138bc8c744e5a453dd7ae22cb7f5608ca9892f19bf98a597becf22eae39216')
        records = [t for p in s['passes'] for t in p['trace']]
        self.assertEqual([t['iteration'] for t in records], list(range(38))+[0])
        self.assertIn('constraint_jacobian', records[0])
        self.assertIn('restoration', records[0])
        self.assertEqual(len(s['binary_cuboids']['cuboids']), 72)
        f = r.artifact.facts
        np.testing.assert_allclose(f['normalized']['moments'],
            [503/512,0,0,-9/4096,32687/393216,32687/393216,32651/393216,0,0,0], atol=1e-10, rtol=0)
        self.assertAlmostEqual(f['physical']['mass'], 503/64)
        self.assertEqual(f['clearance']['squared_distance_exact'], '25/64')
        self.assertEqual(f['newly_carved_components'], 1)
        np.testing.assert_allclose(f['normalized']['central_inertia'],
            np.diag([16432021/98893824,16432021/98893824,32687/196608]),atol=1e-10,rtol=0)
        for a, b in zip(r.outer, self.p.original.arrays()):
            self.assertEqual(a.tobytes(), b.tobytes())
        self.assertEqual(r.identity['geometry_frame'], 'MODEL_COORDINATES_NATIVE_Z_ORIGIN_ZERO')

    def test_all_affine_cross_terms_use_model_units(self):
        transformed = self.m._normalize_model_moments(
            [503/64,503/512,-503/256,5515/512,8549/3072,38723/12288,53513/3072,-5515/2048,5515/4096,-503/2048],
            {'origin':[1/8,-1/4,11/8],'L':2.,'R':1.})
        np.testing.assert_allclose(transformed,
            [503/512,0,0,-9/4096,32687/393216,32687/393216,32651/393216,0,0,0],atol=1e-12,rtol=0)

    def test_repeat_density_and_copy_ownership(self):
        second = self.m.solve_and_seal_source_spin2014(self.p, **self.cfg)
        self.assertEqual(second.scientific_digest, self.solved.scientific_digest)
        double = self.m.solve_and_seal_source_spin2014(problem(2.), **self.cfg)
        np.testing.assert_allclose(double.artifact.facts['normalized']['moments'],
                                   self.solved.artifact.facts['normalized']['moments'], atol=1e-10)
        np.testing.assert_allclose(double.identity['native_run']['geometry']['integrated_moments'],
            2*np.array(self.solved.identity['native_run']['geometry']['integrated_moments']), atol=1e-10)
        self.assertNotEqual(double.scientific_digest, self.solved.scientific_digest)
        self.assertEqual(self.solved.evidence['before'], self.solved.evidence['after'])
        self.assertEqual(self.solved.evidence['scope'],'OBSERVED_FILES_ON_DISK_ONLY')
        copied = self.solved.identity
        copied['native_run']['solve']['status'] = 'FORGED'
        self.assertEqual(self.solved.status, 'BINARY_MOMENT_VALIDATED')
        v, f = self.solved.material; v[:] = 100; f[:] = 0
        self.assertFalse(np.all(self.solved.material[0] == 100))
        with self.assertRaises((FrozenInstanceError, AttributeError)):
            self.solved.status = 'FORGED'
        with self.assertRaises(TypeError): self.m.SourceSpinMaterialResult()

    def test_source_surface_dual_density_exports_both_materials(self):
        dual = SpinArtifactProblem(original=self.p.original.arrays(),
            fixed_regions=[('shell', self.p.fixed_regions[0].mesh.arrays(), 1.25)],
            design_regions=[('design', self.p.design_regions[0].mesh.arrays(), 6.5)],
            length_unit='m', mass_unit='kg', target=self.p.target, minimum_outer_wall=.5)
        result = self.m.solve_and_seal_source_spin2014(dual, **self.cfg)
        self.assertEqual(result.status, 'BINARY_MOMENT_VALIDATED')
        self.assertIsNotNone(result.artifact)
        self.assertIsNotNone(result.boundary_material)
        self.assertIsNotNone(result.interior_material)
        geometry = result.identity['native_run']['geometry']
        self.assertTrue(geometry['material_components_exported'])
        self.assertTrue(geometry['dual_material_partition']['complete_partition'])
        self.assertEqual(geometry['boundary_density'], 1.25)
        self.assertEqual(geometry['interior_density'], 6.5)
        boundary_moments = np.asarray(geometry['boundary_material_moments'])
        interior_moments = np.asarray(geometry['interior_material_moments'])
        integrated_moments = np.asarray(geometry['integrated_moments'])
        # Input shell volume is 7.5 m^3, so its 1.25 density gives 9.375 kg.
        # The retained binary interior is 0.359375 m^3; 6.5 is applied once.
        self.assertAlmostEqual(boundary_moments[0], 9.375, places=12)
        self.assertAlmostEqual(interior_moments[0], 2.3359375, places=12)
        np.testing.assert_allclose(boundary_moments + interior_moments,
                                   integrated_moments, atol=1e-10, rtol=0)
        np.testing.assert_allclose(geometry['integrated_moments'],
            geometry['source_subtraction_moments'], atol=1e-10, rtol=0)
        independent = result.identity['dual_material_independent_audit']
        self.assertEqual(independent['method'], 'INDEPENDENT_SIGNED_TETRAHEDRA_V1')
        self.assertTrue(independent['density_applied_once'])
        self.assertTrue(independent['exact_partition']['complete_partition'])
        self.assertLessEqual(max(independent['normalized_residuals'].values()), 1e-10)
        np.testing.assert_allclose(independent['boundary_moments'], boundary_moments,
                                   atol=1e-10, rtol=0)
        np.testing.assert_allclose(independent['interior_moments'], interior_moments,
                                   atol=1e-10, rtol=0)

        # The independent gate must reject a density applied with the wrong
        # material label, even when the native geometry claim remains intact.
        boundary_snapshot = self.m.MeshSnapshot(result.boundary_material)
        interior_snapshot = self.m.MeshSnapshot(result.interior_material)
        with self.assertRaisesRegex(ValueError, 'dual-material independent mass-moment audit failed'):
            self.m._dual_material_independent_audit(
                boundary_snapshot, interior_snapshot, 2.5, 6.5, geometry, 2.0,
            )

        # A missing phase is not equivalent to a zero-density phase.  It must
        # stop before a combined result can be reported as validated.
        with self.assertRaisesRegex(ValueError, 'cannot integrate an empty material component'):
            self.m._dual_material_independent_audit(
                boundary_snapshot, self.m.MeshSnapshot(None), 1.25, 6.5,
                geometry, 2.0,
            )

    def test_equal_density_source_surface_keeps_single_material_path(self):
        # The equal-density limit must use the ordinary single-material path;
        # splitting a region must not multiply or otherwise alter its density.
        result = self.solved
        self.assertIsNone(result.boundary_material)
        self.assertIsNone(result.interior_material)
        self.assertIsNone(result.identity['dual_material_independent_audit'])
        geometry = result.identity['native_run']['geometry']
        np.testing.assert_allclose(geometry['integrated_moments'],
                                   geometry['source_subtraction_moments'],
                                   atol=1e-10, rtol=0)

    def test_resource_and_wall_failures_keep_numerical_status(self):
        cfg = {**self.cfg, 'boundary_settings': self.m.SourceBoundarySettings(maximum_intermediate_complexity=4)}
        limited = self.m.solve_and_seal_source_spin2014(self.p, **cfg)
        self.assertEqual(limited.status, 'BINARY_MOMENT_VALIDATED')
        self.assertEqual(limited.identity['native_run']['geometry']['status'], 'RESOURCE_LIMIT_UNKNOWN')
        self.assertIsNone(limited.artifact)
        wall = self.m.solve_and_seal_source_spin2014(problem(wall=np.nextafter(5/8, 1.)), **self.cfg)
        self.assertEqual(wall.status, 'BINARY_MOMENT_VALIDATED')
        self.assertIsNone(wall.artifact)

    def test_quarter_fill_and_all_empty_objectives_do_not_seal(self):
        quarter = problem(outer=box((-.5+1/8192,-.5,0),(.5+1/8192,.5,1)),
                          design=box((.125,-.125,.625),(.375,.125,.875)))
        cfg = {**self.cfg, 'source_settings': replace(self.cfg['source_settings'], maximum_depth=5),
               'solver_settings': replace(self.cfg['solver_settings'], maximum_optimization_level=5)}
        r = self.m.solve_and_seal_source_spin2014(quarter, **cfg)
        self.assertEqual(r.status, 'BINARY_VALIDATION_FAILED'); self.assertIsNone(r.artifact)
        for objective, mode in [('STATIC_LOWERING','STATIC_BALANCE'),('YOYO','YO_YO'),
            ('TOP_INERTIA_ONLY','SPINNING_TOP'),('TOP_LOWERING_ONLY','SPINNING_TOP'),('TOP_MASS_AND_LOWERING','SPINNING_TOP')]:
            p = problem(outer=box((-.5,)*3,(.5,)*3), design=box((.39,)*3,(.40,)*3), mode=mode)
            cfg = {**self.cfg, 'source_settings': self.m.SourceGeometrySettings(root_origin=(-1,-1,-1),root_side=2,maximum_depth=3),
                'solver_settings': replace(self.cfg['solver_settings'], maximum_optimization_level=3,
                    objective=objective, initial_phi=0 if mode=='STATIC_BALANCE' else .37)}
            r = self.m.solve_and_seal_source_spin2014(p, **cfg)
            self.assertEqual(r.status,'RELAXATION_UNKNOWN'); self.assertIsNone(r.artifact)
            self.assertIsNone(r.identity['native_run']['solve']['passes'][0]['final_kkt']['full'])

    def test_settings_and_target_rejections(self):
        for cls, key, bad in [(self.m.SourceGeometrySettings,'maximum_depth',True),
            (self.m.SourceGeometrySettings,'root_side',2**53+1),
            (self.m.SourceGeometrySettings,'root_origin',(False,0,0)),
            (self.m.SourceSolverSettings,'initial_fill',float('nan')),
            (self.m.SourceSolverSettings,'objective','BAD'),
            (self.m.SourceSolverSettings,'formulation','LEGACY'),
            (self.m.AuditedSlqpSettings,'algorithm','LEGACY_V3'),
            (self.m.AuditedSlqpSettings,'lp_iteration_limit',False),
            (self.m.SourceBoundarySettings,'maximum_output_vertices',0)]:
            with self.subTest(key=key), self.assertRaises((ValueError,TypeError)): cls(**{key:bad})
        with self.assertRaises(TypeError): self.m.SourceBoundarySettings(minimum_outer_wall=.5)
        for t in ({'axis':[1,0,0]}, {'origin':[0,0,1]}, {'support_points':[[0,0,0],[0,0,0]]}):
            with self.assertRaisesRegex(ValueError,'NOT_SUPPORTED'):
                self.m.solve_and_seal_source_spin2014(problem(**t), **self.cfg)
        with self.assertRaises(ValueError): problem(axis=[0,0,0])
        with self.assertRaisesRegex(ValueError,'NOT_SUPPORTED'):
            self.m.solve_and_seal_source_spin2014(problem(mode='BUOYANT_STABILITY'), **self.cfg)
        o = self.p.original
        fixed = [(r.region_id,r.mesh,r.density) for r in self.p.fixed_regions]
        for design in [[('left',box((-.5,-.5,1),(0,.5,1.5)),1.),
                        ('right',box((0,-.5,1),(.5,.5,1.5)),1.)]]:
            p = SpinArtifactProblem(original=o,fixed_regions=fixed,design_regions=design,
                length_unit='m',mass_unit='kg',target=self.p.target,minimum_outer_wall=.5)
            with self.assertRaisesRegex(ValueError,'NOT_SUPPORTED'):
                self.m.solve_and_seal_source_spin2014(p,**self.cfg)
        with self.assertRaises(ValueError):
            SpinArtifactProblem(original=o,fixed_regions=fixed,
                design_regions=[('design',self.p.design_regions[0].mesh,1.)],
                length_unit='m',mass_unit='kg',target=self.p.target,minimum_outer_wall=.5,
                validation_scope='DEFORMATION')

    def test_source_target_frame_binds_problem_axis_and_native_settings(self):
        axes=((0.,0.,1.),(0.,1.,0.),(-1.,0.,0.))
        solver=replace(self.cfg['solver_settings'], target_axes=axes)
        rotated_target=problem(axis=[1,0,0])
        result=self.m.solve_and_seal_source_spin2014(
            rotated_target, source_settings=self.cfg['source_settings'],
            solver_settings=solver, boundary_settings=self.cfg['boundary_settings'],
        )
        settings=result.identity['native_run']['settings']['solver_settings']
        self.assertEqual(settings['target_axes'],[list(row) for row in axes])
        self.assertEqual(settings['target_origin'],[0.,0.,0.])
        with self.assertRaisesRegex(ValueError,'NOT_SUPPORTED'):
            self.m.solve_and_seal_source_spin2014(
                problem(axis=[0,0,1]), source_settings=self.cfg['source_settings'],
                solver_settings=solver, boundary_settings=self.cfg['boundary_settings'],
            )

    def test_direct_native_lossless_mesh_and_settings_boundary(self):
        o, d = self.p.original.arrays(), self.p.design_regions[0].mesh.arrays()
        args = dict(source_settings=asdict(self.cfg['source_settings']), solver_settings=asdict(self.cfg['solver_settings']),
            boundary_settings=asdict(self.cfg['boundary_settings']), density=1., minimum_outer_wall=.5)
        badv = o[0].tolist(); badv[0][0] = False
        lossy = o[0].tolist(); lossy[0][0] = 2**53+1
        badf = o[1].tolist(); badf[0][0] = True
        for v,f in [(badv,o[1]),(lossy,o[1]),(o[0],badf),(o[0],o[1].astype(float))]:
            with self.assertRaises((ValueError,TypeError)):
                self.n.solve_source_spin2014_material(v,f,*d,**args)
        for field, bad in [('maximum_depth',True),('root_side',2**53+1),('extra',1)]:
            a = {**args, 'source_settings': {**args['source_settings'],field:bad}}
            with self.assertRaises((ValueError,TypeError)):
                self.n.solve_source_spin2014_material(*o,*d,**a)
        for key,bad in [('density',True),('minimum_outer_wall',np.longdouble(.5)),('density',float('inf'))]:
            with self.assertRaises((ValueError,TypeError)):
                self.n.solve_source_spin2014_material(*o,*d,**{**args,key:bad})
        badslqp={**args['solver_settings']['slqp'],'lp_iteration_limit':True}
        with self.assertRaises((ValueError,TypeError)):
            self.n.solve_source_spin2014_material(*o,*d,**{**args,'solver_settings':{**args['solver_settings'],'slqp':badslqp}})
        with self.assertRaises(TypeError): self.n.NativeSourceSpinMaterialResult()
        direct=self.n.solve_source_spin2014_material(*o,*d,**args)
        identity=direct.diagnostics
        storage=direct.trace_storage
        compact_lines=list(direct.iter_stage_audit_json())
        compact_events=[__import__('json').loads(line) for line in compact_lines]
        self.assertEqual(compact_events[0]['kind'],'header')
        self.assertEqual(compact_events[0]['trace_counts'],
                         [len(p['trace']) for p in identity['solve']['passes']])
        iteration_events=[event for event in compact_events if event['kind']=='iteration']
        self.assertEqual(len(iteration_events),storage['record_count'])
        counts=compact_events[0]['trace_counts']; cursor=0
        for pass_index,count in enumerate(counts):
            pass_events=[event for event in iteration_events if event['pass_index']==pass_index]
            self.assertEqual([event['trace_index'] for event in pass_events],list(range(count)))
            for event in pass_events[:-1]: self.assertEqual(event['payload'],{})
            self.assertTrue(pass_events[-1]['payload'].get('final_state'))
            self.assertIn('state',pass_events[-1]['payload'])
            cursor+=count
        self.assertEqual(cursor,storage['record_count'])
        self.assertEqual(identity['settings']['source_settings'], {**args['source_settings'],'root_origin':[-1.,-1.,0.]})
        self.assertEqual(
            identity['settings']['solver_settings'],
            json.loads(json.dumps(args['solver_settings'])),
        )
        self.assertEqual(identity['settings']['boundary_settings'],args['boundary_settings'])
        direct.outer[0][:]=100
        direct.material[1][:]=0
        direct.diagnostics['solve']['status']='FORGED'
        self.assertEqual(direct.diagnostics,identity)
        self.assertEqual(direct.outer[0].tobytes(),o[0].tobytes())

    def test_persisted_bundle_reintegrates_and_rejects_tampering(self):
        from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
        with tempfile.TemporaryDirectory(prefix='task-5f-bundle-') as tmp:
            p = Path(tmp)/'sealed'
            saved = save_spin_artifact(self.solved.artifact, p)
            self.assertEqual(saved.scientific_digest, self.solved.artifact.scientific_digest)
            with self.assertRaises(ValueError): load_spin_artifact(p, expected_problem_digest='0'*64)
            mesh = next(p.glob('mesh_*_vertices.npy'))
            data = bytearray(mesh.read_bytes()); data[-1] ^= 1; mesh.write_bytes(data)
            with self.assertRaises(ValueError): load_spin_artifact(p, expected_problem_digest=self.p.input_digest)


class SourceRunnerTests(unittest.TestCase):
    def test_transition_runner_preserves_real_split_and_binary_failure(self):
        root = Path(__file__).resolve().parents[2]
        path = root/'experiments/spin_it_reproduction/strict/run_source_spin2014_adaptive_transitions.py'
        self.assertTrue(path.is_file(), 'missing source-derived adaptive transition runner')
        spec = importlib.util.spec_from_file_location('source_transition_runner',path)
        runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
        import pybind11
        import _make_it_stand as native
        from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        tmp = Path(tempfile.mkdtemp(prefix='source-transition-runner-test-',dir=root/'build'))
        summary = runner.run(tmp/'transition',
            build_directory=Path(native.__file__).resolve().parents[1],
            dependency_roots={'pybind11':Path(pybind11.get_include())})
        self.assertEqual(summary['status'], 'BINARY_VALIDATION_FAILED')
        self.assertEqual(summary['classification'], 'SYNTHETIC_REGRESSION')
        self.assertFalse(summary['sealed'])
        native_run = __import__('json').loads((tmp/'transition'/'native-run.json').read_text())
        solve = native_run['solve']
        self.assertEqual(solve['termination_policy'], 'PAPER_POST_TRANSITION_V1')
        self.assertEqual(len(solve['passes']), 3)
        self.assertTrue(all(item['status'] == 'CONVERGED' for item in solve['passes']))
        self.assertGreater(sum(len(item['split']['changes']) for item in solve['passes']), 0)
        self.assertEqual(sum(len(item['merge']['changes']) for item in solve['passes']), 0)
        self.assertGreater(native_run['source']['stats']['boundary_nodes'], 0)
        self.assertTrue(solve['passes'][-1]['paper_resolution_complete'])
        inspection = __import__('json').loads((tmp/'transition'/'domain-inspection.json').read_text())
        self.assertEqual(inspection['source']['source_hash'], native_run['source']['source_hash'])
        self.assertEqual(inspection['source']['partition_hash'], native_run['source']['partition_hash'])
        self.assertGreater(len(inspection['source_cells']), 0)
        self.assertGreater(len(inspection['edges']), 0)
        evidence = PaperRunEvidenceV2.from_json((tmp/'transition'/'evidence.json').read_text())
        self.assertEqual(check_freshness(evidence).status, 'CURRENT')

    def test_runner_persists_actual_success_failure_and_freshness(self):
        root = Path(__file__).resolve().parents[2]
        path = root/'experiments/spin_it_reproduction/strict/run_source_spin2014_material.py'
        self.assertTrue(path.is_file(), 'missing actual synthetic source runner')
        spec = importlib.util.spec_from_file_location('source_runner',path)
        runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
        import pybind11
        import _make_it_stand as native
        from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        from make_it_stand.verification.spin_bundle import load_spin_artifact
        from make_it_stand.verification.source_spin_material import SourceBoundarySettings
        # Evidence remains on disk, including failed attempts, for review.
        tmp = Path(tempfile.mkdtemp(prefix='task-5f-runner-test-',dir=root/'build'))
        kwargs = {'build_directory':Path(native.__file__).resolve().parents[1],
                  'dependency_roots':{'pybind11':Path(pybind11.get_include())}}
        summary = runner.run(tmp/'success', **kwargs)
        self.assertEqual(summary['status'],'BINARY_MOMENT_VALIDATED')
        self.assertTrue(summary['sealed'])
        payload = (tmp/'success'/'scientific.json').read_text()
        data = __import__('json').loads(payload)
        artifact = load_spin_artifact(tmp/'success'/'sealed',expected_problem_digest=data['identity']['problem_digest'])
        self.assertEqual(artifact.scientific_digest,data['identity']['artifact_scientific_digest'])
        evidence = PaperRunEvidenceV2.from_json((tmp/'success'/'evidence.json').read_text())
        self.assertEqual(check_freshness(evidence).status,'CURRENT')
        with self.assertRaises(FileExistsError): runner.run(tmp/'success', **kwargs)
        failed = runner.run(tmp/'failure',boundary_settings=SourceBoundarySettings(maximum_intermediate_complexity=4),**kwargs)
        self.assertFalse(failed['sealed'])
        bad = PaperRunEvidenceV2.from_json((tmp/'failure'/'evidence.json').read_text()).as_dict()
        self.assertEqual(bad['scientific']['events'][-1]['status'],'failed')
        self.assertFalse((tmp/'failure'/'sealed').exists())
        with (tmp/'success'/'inputs.json').open('a') as stream: stream.write(' ')
        self.assertEqual(check_freshness(evidence).status,'STALE_RESULT')


if __name__ == '__main__':
    unittest.main()
