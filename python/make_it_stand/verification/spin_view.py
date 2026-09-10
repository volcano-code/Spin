"""Copy-owned model-frame scenes admitted only through the strict bundle loader.

This adapter does not run a solver, alter material, or admit solver execution.
Eigenspaces describe complete density-weighted world central inertia; axis signs
are conventional. Near-repeated eigenvalues are conservatively grouped at a
relative 1e-10 display tolerance, not asserted to have distinct physical axes.
"""
import numpy as np

from .spin_bundle import load_spin_artifact


def _principal_axes(inertia):
    values, vectors = np.linalg.eigh(np.asarray(inertia, dtype=np.float64))
    tolerance = 1e-10 * float(np.max(np.abs(values)))
    groups = []
    for i, value in enumerate(values):
        if groups and value - values[groups[-1][0]] <= tolerance:
            groups[-1].append(i)
        else:
            groups.append([i])
    axes = []
    for group in groups:
        if len(group) == 1:
            i = group[0]
            direction = vectors[:, i].copy()
            if direction[np.argmax(np.abs(direction))] < 0:
                direction *= -1
            axes.append({'eigenvalue': float(values[i]), 'direction': direction.tolist()})
    return {'eigenvalues': values.tolist(), 'multiplicities': [len(g) for g in groups],
            'unique_axes': axes, 'relative_display_tolerance': 1e-10,
            'sign_convention': 'LARGEST_COMPONENT_POSITIVE'}


def load_spin_scene(bundle_directory, *, expected_problem_digest):
    """Independently admit a bundle, then return owned arrays and computed facts."""
    artifact = load_spin_artifact(bundle_directory, expected_problem_digest=expected_problem_digest)
    p, facts = artifact.problem, artifact.facts
    meshes = []

    def add(key, kind, mesh, density=None):
        vertices, faces = mesh.arrays()
        meshes.append({'key': key, 'kind': kind, 'vertices': vertices, 'faces': faces,
                       'density': density, 'identity': mesh.identity})

    add('original', 'original', p.original)
    add('outer', 'outer', artifact.outer)
    for r in p.fixed_regions:
        add('fixed:' + r.region_id, 'fixed', r.mesh, r.density)
    add('material', 'material', artifact.final_material)
    density = {r.region_id: r.density for r in p.design_regions}
    for region_id, retained, _ in artifact.design_results:
        add('retained:' + region_id, 'retained', retained, density[region_id])
    for r in p.design_regions:
        add('design:' + r.region_id, 'design', r.mesh, r.density)
    add('cavity', 'cavity', artifact.newly_carved)
    disclosures = {name: facts[name] for name in
                   ('project_status', 'original_evidence_admission', 'physical_validation',
                    'build_source_attestation', 'solver_target_feasibility')}
    disclosures['numerical_solver_trace_convergence'] = 'NOT_AUTHENTICATED_BY_GEOMETRY_BUNDLE'
    return {'meshes': meshes, 'physical': facts['physical'],
            'principal_axes': _principal_axes(facts['physical']['central_inertia']),
            'target': p.target, 'support_display': 'TARGET_POINTS_ONLY_NO_CERTIFIED_HULL',
            'geometry_status': facts['geometry_status'], 'artifact_digest': artifact.scientific_digest,
            'evidence_digest': artifact.evidence_digest, 'problem_digest': p.input_digest,
            'validation_scope': facts['validation_scope'], 'disclosures': disclosures,
            'length_unit': p.length_unit, 'mass_unit': p.mass_unit}
