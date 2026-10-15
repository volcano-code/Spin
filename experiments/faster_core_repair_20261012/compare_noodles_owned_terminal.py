"""Actual saved noodles terminal: independent exact representation comparison.

No optimization, source reconstruction or input-partition re-preparation.
This diagnostic cannot grant the full material/case seal.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import trimesh

from make_it_stand import spin
from make_it_stand.verification.cherchi_run_evidence import load_cherchi_terminal_prefix
from make_it_stand.verification.cherchi_run_journal import RunJournal
from make_it_stand.verification.spin_problem import MeshSnapshot
from make_it_stand.verification import faster_material as fm

HERE = Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    load = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    saved = HERE/'noodles-terminal-recovery-job-01/output'
    frozen_path = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    frozen = load(frozen_path)
    original = HERE/'noodles-formal-cherchi-seed-0-job-01'
    with RunJournal(output, durability_profile='PER_EVENT_FSYNC_V1', json_value=fm._json_value) as journal:
        def stage(name, operation):
            start = time.monotonic_ns()
            journal.emit(dict(kind='STAGE_STARTED', stage=name, monotonic_ns=start))
            print('STAGE_STARTED', name, flush=True)
            try:
                result = operation()
            except BaseException as error:
                journal.emit(dict(kind='STAGE_FAILED', stage=name, error=f'{type(error).__name__}: {error}',
                                  elapsed_ns=time.monotonic_ns()-start))
                raise
            elapsed = (time.monotonic_ns()-start)/1e9
            journal.emit(dict(kind='STAGE_RETURNED', stage=name, elapsed_seconds=elapsed))
            print('STAGE_RETURNED', name, elapsed, flush=True)
            return result, elapsed
        prefix, _ = stage('verify_original_accepted_prefix', lambda:
                          load_cherchi_terminal_prefix(original/'output/run'))
        original_registration = load(original/'output/registration.json')
        if sha(frozen_path) != original_registration['request_sha256']:
            raise ValueError('FROZEN_REQUEST_CHANGED')
        for name, expected in frozen['input_sources'].items():
            if sha(Path(name)) != expected:
                raise ValueError('FROZEN_INPUT_CHANGED:'+name)
        recovery_registration = load(saved/'registration.json')
        if recovery_registration['accepted_material'] != prefix['terminal_reference']:
            raise ValueError('RECOVERED_TERMINAL_REFERENCE_CHANGED')
        carved = load(saved/'terminal-complement.json')
        if carved['quadric_source_identity'] != prefix['material']['quadric_source_identity']:
            raise ValueError('RECOVERED_COMPLEMENT_SOURCE_MISMATCH')
        request_data = load(saved/'prepared-request.json')
        request_data['settings'] = spin.SpinSolverConfig(**request_data['settings'])
        request = spin.SpinItFasterRequest(**request_data)
        normalization = frozen['normalization']
        length, density = float(normalization['L']), float(normalization['R'])
        if length != request.original_units_per_normalized_unit or request.material_density != 1.:
            raise ValueError('UNSUPPORTED_COORDINATE_OR_DENSITY_CONTRACT')
        paths = {Path(name).name: Path(name) for name in frozen['input_sources']}
        mesh = trimesh.load(paths['outer-si-candidate.obj'], force='mesh', process=False)
        actual_origin = mesh.vertices.min(0)+(mesh.vertices.max(0)-mesh.vertices.min(0))*.5
        if (actual_origin.tolist() != normalization['origin'] or
                float(max(mesh.vertices.max(0)-mesh.vertices.min(0))) != length):
            raise ValueError('ORIGINAL_NORMALIZATION_CHANGED')
        journal.publish('registration.json', dict(schema='noodles-owned-terminal-comparison-v1',
            original_job=str(original), original_run_complete=False, case_complete=False,
            accepted_reference=prefix['terminal_reference'], complement_sha256=sha(saved/'terminal-complement.json'),
            frozen_sha256=sha(frozen_path), native_sha256=sha(Path(spin._native_module().__file__)),
            optimization_invocations=0, source_reconstruction_invocations=0))
        outer, outer_seconds = stage('original_exact_outer', lambda: MeshSnapshot((mesh.vertices, mesh.faces)))
        normalized = fm._target_rationals_to_normalized_model(carved['variable_rational_vertices'], request)
        world = [[fm.fraction_to_rational_text(fm.fraction_from_rational_text(value)*Fraction(length)+
                   Fraction(float(normalization['origin'][axis]))) for axis, value in enumerate(row)]
                 for row in normalized]
        removed, carved_seconds = stage('original_exact_carved', lambda:
            MeshSnapshot.from_rationals(world, carved['variable_faces']))
        ov, of = outer.arrays()
        cv, cf = removed.arrays()
        exact_outer = [[fm.fraction_to_rational_text(Fraction(float(value))) for value in row] for row in ov]
        combined = exact_outer+world
        faces = np.vstack((of, cf[:, ::-1]+len(ov)))
        baseline, baseline_seconds = stage('baseline_complete_boundary_reconstruction', lambda:
            MeshSnapshot.from_rationals(combined, faces))
        journal.publish('baseline-facts.json', baseline.exact_region().facts)
        region, difference_seconds = stage('owned_exact_interior_difference', lambda:
            spin._native_module().make_exact_interior_difference_region(outer.exact_region(), removed.exact_region()))
        if not region.facts['completed'] or not region.facts['valid']:
            journal.publish('candidate-failure.json', region.facts)
            raise ValueError('OWNED_DIFFERENCE_NOT_VALID')
        candidate, export_seconds = stage('owned_boundary_export', lambda: MeshSnapshot.from_exact_region(region))
        equal_boundary = (candidate.rational_vertices == baseline.rational_vertices and
                          all(a.tobytes() == b.tobytes() for a, b in zip(candidate.arrays(), baseline.arrays())))
        if not equal_boundary:
            raise ValueError('OWNED_EXPORT_CHANGED_ACTUAL_BOUNDARY')
        relation, relation_seconds = stage('independent_exact_set_comparison', lambda:
            spin._native_module().relate_exact_material_regions(baseline.exact_region(), region))
        journal.publish('exact-set-relation.json', relation)
        if not relation['completed'] or not relation['equal']:
            raise ValueError('OWNED_EXPORT_CHANGED_ACTUAL_REGION')
        moments, integration_seconds = stage('independent_world_material_integration', lambda:
                                             fm._raw_moments(candidate, density))
        retained = load(saved/'independent-material-moments.json')['retained']
        expected = fm._transform_moments(retained, request.target_axes, request.target_origin)
        expected = np.asarray(fm._transform_moments(expected, np.eye(3)*length, normalization['origin']))
        expected *= density*length**3
        scales = density*np.asarray([length**power for power in (3,4,4,4,5,5,5,5,5,5)])
        error = float(np.max(np.abs(np.asarray(moments)-expected)/scales))
        if error > 1e-6:
            raise ValueError('SAVED_TERMINAL_TOTAL_MOMENTS_MISMATCH')
        journal.publish('report.json', dict(passed=True, equal_boundary=True, exact_relation=relation,
            normalized_total_moment_error=error, outer_seconds=outer_seconds, carved_seconds=carved_seconds,
            baseline_reconstruction_seconds=baseline_seconds, candidate_difference_seconds=difference_seconds,
            candidate_export_seconds=export_seconds, relation_seconds=relation_seconds,
            independent_integration_seconds=integration_seconds, optimization_invocations=0,
            full_material_seal=False, case_complete=False,
            scope='Exact same saved terminal representation comparison; not complete partition/seal or a new optimization'))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
