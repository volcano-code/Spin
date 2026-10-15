"""Exact fixed-shell and interior-domain check for pure hollowing.

Spin-It (2014) Section 5.1 / Figure 3 restricts removal to an interior
domain and retains its boundary shell. This sidecar checks that construction
against a locked input, including clearance of the WHOLE design domain, not
only the selected cavity. It does not measure every original thin feature,
integrate mass, admit a material to physics, or authenticate author assets.

No existing sealed-artifact format or historical facts are changed.
"""
from .spin_common import native, require
from .spin_problem import MeshSnapshot, SpinArtifactProblem, partition


def check_fixed_boundary_shell(problem, *, outer, design_results,
                               final_material, newly_carved):
    require(isinstance(problem, SpinArtifactProblem), 'locked SpinArtifactProblem required')
    require(bool(problem.fixed_regions), 'explicit fixed boundary shell required')
    outer = MeshSnapshot(outer)
    material, carved = MeshSnapshot(final_material), MeshSnapshot(newly_carved)
    require(outer == problem.original, 'outer canonical geometry bytes changed')
    require(not material.explicit_empty, 'final material must be nonempty')
    original_region = outer.exact_region()
    fixed = [r.mesh.exact_region() for r in problem.fixed_regions]
    interiors = {r.region_id: r.mesh.exact_region() for r in problem.design_regions}
    initial = partition(original_region, fixed + list(interiors.values()))

    clearances = []
    for region_id, interior in interiors.items():
        clearance = native().check_exact_boundary_clearance(
            original_region, interior, problem.minimum_outer_wall)
        require(clearance['completed'] and clearance['applicable'] and
                clearance['at_least_minimum'], 'design-domain clearance violates fixed shell')
        clearances.append({'id': region_id, 'clearance': clearance})

    require(isinstance(design_results, (list, tuple)), 'explicit design-result list required')
    seen, results, kept_regions, removed_regions, identities = set(), [], [], [], []
    for row in design_results:
        require(isinstance(row, (tuple, list)) and len(row) == 3,
                'design result must contain id, retained mesh, carved mesh')
        region_id, kept, removed = row
        require(isinstance(region_id, str) and region_id in interiors and region_id not in seen,
                'duplicate or unknown design ID')
        seen.add(region_id)
        kept, removed = MeshSnapshot(kept), MeshSnapshot(removed)
        kr, rr = kept.exact_region(), removed.exact_region()
        results.append({'id': region_id, 'partition': partition(interiors[region_id], [kr, rr])})
        kept_regions.append(kr)
        removed_regions.append(rr)
        identities.append({'id': region_id, 'retained': kept.identity, 'removed': removed.identity})
    require(seen == set(interiors), 'missing design ID')
    final_partition = partition(material.exact_region(), fixed + kept_regions)
    carved_partition = partition(carved.exact_region(), removed_regions)
    return {
        'schema': 'spin2014-fixed-boundary-shell-check-v1',
        'status': 'FIXED_BOUNDARY_SHELL_VALIDATED',
        'scope': 'HOLLOW_ONLY_LOCKED_INPUT_DOMAIN_AND_PARTITIONS',
        'problem_digest': problem.input_digest,
        'minimum_domain_clearance': problem.minimum_outer_wall,
        'length_unit': problem.length_unit,
        'domain_clearances': clearances,
        'initial_partition': initial, 'design_partitions': results,
        'final_material_partition': final_partition, 'carved_partition': carved_partition,
        'output_identity': {'outer': outer.identity, 'material': material.identity,
                            'carved': carved.identity, 'design_results': identities},
        'full_feature_thickness': 'UNKNOWN', 'material_admission': 'NOT_GRANTED',
        'moment_validation': 'NOT_PERFORMED',
        'original_paper_reproduction': 'NOT_GRANTED',
    }
