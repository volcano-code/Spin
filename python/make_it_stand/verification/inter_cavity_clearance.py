"""Exact separation of distinct bounded voids in an authenticated material.

This is a supplemental observation, not an extension of the frozen spin seal.
The caller supplies an explicit separation threshold; minimum_outer_wall is
not silently reinterpreted. Complete pairwise boundary separation is only one
part of a full wall/feature audit. A short segment can cross exterior air, so a
below-threshold observation is NOT by itself a material-thickness failure.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

from .spin_common import digest, implementation_observation, native, require, scalar
from .spin_problem import MeshSnapshot
from .run_evidence import _record


def audit_inter_cavity_clearance(mesh, *, minimum_clearance, maximum_pairs):
    """Return recomputed bounded-void distances using exact native predicates.

    ``mesh`` is a MeshSnapshot, a (vertices, faces) pair, or explicit None.
    A void boundary consists of its negative shell AND its immediate positive
    child shells (material islands). Face/vertex IDs in witnesses refer to the
    original input, despite reversal/compaction of temporary void copies.

    ``maximum_pairs`` bounds pair enumeration, not total native geometry work;
    large inputs still require the project's supervised resource-limited job.
    Exhausted/incomplete checks publish no minimum or successful comparison.
    """
    threshold = scalar(minimum_clearance, 'minimum inter-cavity clearance')
    require(threshold >= 0, 'minimum inter-cavity clearance must be nonnegative')
    require(type(maximum_pairs) is int and 0 < maximum_pairs <= 2**31 - 1,
            'maximum_pairs must be a positive signed32 integer')
    snapshot = mesh if isinstance(mesh, MeshSnapshot) else MeshSnapshot(mesh)
    n = native()
    material = snapshot.exact_region()
    facts = material.facts
    before = implementation_observation()
    before['files']['inter_cavity_clearance'] = _record(Path(__file__), label='inter_cavity_clearance')
    count = facts['bounded_void_components']
    required = count * (count - 1) // 2
    result = {
        'schema': 'inter-cavity-clearance-v1',
        'scope': 'ALL_DISTINCT_BOUNDED_VOID_BOUNDARIES_ONLY',
        'input': snapshot.identity,
        'material_region_digest': facts['region_digest'],
        'minimum_clearance_exact': str(Fraction(threshold)),
        'maximum_pairs': maximum_pairs,
        'cavity_count': count,
        'pairs_required': required,
        'pairs_checked': 0,
        'cavities': [],
        'completed': False,
        'state': 'INCOMPLETE',
        'error': '',
        'at_least_minimum': None,
        'squared_distance_exact': None,
        'witness_cavities': None,
        'witness_faces': None,
        'complete_wall_thickness': 'UNKNOWN',
        'manufacturing_admission': 'NOT_GRANTED',
        'physical_admission': 'NOT_GRANTED',
        'implementation': before,
    }

    def finish():
        after = implementation_observation()
        after['files']['inter_cavity_clearance'] = _record(Path(__file__), label='inter_cavity_clearance')
        require(after == before, 'checker files changed during inter-cavity audit')
        result['audit_digest'] = digest(result)
        return result

    if count < 2:
        result.update(completed=True, state='NOT_APPLICABLE')
        return finish()
    if required > maximum_pairs:
        result.update(state='RESOURCE_LIMIT', error='all cavity pairs exceed maximum_pairs')
        return finish()

    try:
        vertices, faces = snapshot.arrays()
        boundary = n.audit_exact_triangle_boundary(vertices, faces, 'AnyNonzero')
        require(boundary['boundary_valid'] and boundary['input_digest'] == facts['input_digest'],
                'material and boundary identities disagree')
        signs, parents = facts['original_signs'], facts['parents']
        require(signs == boundary['component_volume_signs'], 'shell orientation identities disagree')
        component_ids = np.asarray(boundary['face_components'], dtype=np.int64)
        cavity_ids = [i for i, sign in enumerate(signs) if sign == -1]
        require(len(cavity_ids) == count, 'bounded void count disagrees with shell hierarchy')
        regions, source_faces = [], []
        for cavity in cavity_ids:
            shells = [cavity] + [i for i, parent in enumerate(parents) if parent == cavity]
            require(all(signs[i] == 1 for i in shells[1:]), 'invalid material-island orientation')
            ids = np.flatnonzero(np.isin(component_ids, shells))
            used, compact = np.unique(faces[ids], return_inverse=True)
            # Original mesh is never edited: negative material-side walls become
            # outward void boundaries; positive island shells become inward.
            void_faces = compact.reshape((-1, 3))[:, ::-1].copy()
            region = n.make_exact_material_region(vertices[used], void_faces)
            require(region.facts['completed'] and region.facts['valid'],
                    'exact void reconstruction failed: ' + region.facts['error'])
            require(region.facts['material_components'] == 1,
                    'reconstructed void is not one connected region')
            regions.append(region)
            source_faces.append(ids)
            result['cavities'].append({
                'shell_id': cavity, 'boundary_shells': shells,
                'source_face_ids': ids.tolist(), 'source_vertex_ids': used.tolist(),
                'region_digest': region.facts['region_digest'],
            })

        best = None
        for a, b in combinations(range(count), 2):
            measured = n.check_exact_boundary_clearance(regions[a], regions[b], threshold)
            require(measured['completed'] and measured['applicable'],
                    'exact cavity distance incomplete: ' + measured['error'])
            ia, ib = measured['witness_faces']
            candidate = (Fraction(measured['squared_distance_exact']),
                         cavity_ids[a], cavity_ids[b],
                         int(source_faces[a][ia]), int(source_faces[b][ib]))
            if best is None or candidate < best:
                best = candidate
            result['pairs_checked'] += 1
        require(result['pairs_checked'] == required and best is not None,
                'incomplete cavity pair enumeration')
        squared, a, b, fa, fb = best
        result.update(completed=True, state='COMPLETE',
                      squared_distance_exact=f'{squared.numerator}/{squared.denominator}',
                      at_least_minimum=squared >= Fraction(threshold)**2,
                      witness_cavities=[a, b], witness_faces=[fa, fb])
    except (ValueError, RuntimeError, MemoryError) as exc:
        result.update(state='RESOURCE_LIMIT' if isinstance(exc, MemoryError) else 'INCOMPLETE',
                      error=str(exc))
    return finish()
