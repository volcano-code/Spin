"""Recomputed, native-backed hollow mesh/material validation, not solver admission."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import weakref

import numpy as np

from .mass_moments import TriangleMesh, compute_mass_properties
from .spin_common import (FRAME, ORDER, SCOPE, TOLERANCE, canonical, digest,
                          implementation_observation, keys, native, parse, require, scalar)
from .spin_problem import MeshSnapshot, SpinArtifactProblem, partition
from .rational_text import fraction_from_rational_text, fraction_to_rational_text


def _ten(mass, first, second):
    return [mass, *first, second[0][0], second[1][1], second[2][2],
            second[1][2], second[0][2], second[0][1]]


def _integrate(mesh, density, frame=None):
    rational_vertices = getattr(mesh, 'rational_vertices', None)
    if rational_vertices is not None:
        from .mass_moments import integrate_exact_rational_raw_moments
        vertices = [[fraction_from_rational_text(value) for value in row] for row in rational_vertices]
        faces = mesh.arrays()[1].tolist()
        weight = Fraction(float(density))
        if frame is not None:
            origin = [Fraction(float(value)) for value in frame['origin']]
            length = Fraction(float(frame['L']))
            vertices = [[(value - origin[axis]) / length for axis, value in enumerate(row)]
                        for row in vertices]
            weight /= Fraction(float(frame['R']))
        exact = integrate_exact_rational_raw_moments(vertices, faces)
        weighted = tuple(value * weight for value in exact)
        mass, first = weighted[0], weighted[1:4]
        second = ((weighted[4], weighted[9], weighted[8]),
                  (weighted[9], weighted[5], weighted[7]),
                  (weighted[8], weighted[7], weighted[6]))
        center = tuple(first[axis] / mass for axis in range(3))
        centered = tuple(tuple(second[i][j] - first[i] * first[j] / mass
                                for j in range(3)) for i in range(3))
        trace = sum(centered[i][i] for i in range(3))
        inertia = tuple(tuple((trace if i == j else 0) - centered[i][j]
                              for j in range(3)) for i in range(3))
        geometric = {'volume': float(exact[0]), 'first_moment': [float(x) for x in exact[1:4]],
                     'second_moment': [[float(x) for x in row] for row in
                         ((exact[4], exact[9], exact[8]), (exact[9], exact[5], exact[7]),
                          (exact[8], exact[7], exact[6]))]}
        properties = {'mass': float(mass), 'moments': [float(x) for x in weighted],
                      'center_of_mass': [float(x) for x in center],
                      'central_inertia': [[float(x) for x in row] for row in inertia],
                      'exact_rational_moments': [fraction_to_rational_text(x) for x in weighted]}
        require(properties['mass'] > 0, 'density-weighted region mass underflow')
        canonical(properties)
        return geometric, properties
    v, f = mesh.arrays()
    if frame is not None:
        # Preserve input geometry through the independent affine integration;
        # a rounded vertex transform can collapse genuine Boolean slivers.
        v = np.array([[(Fraction(float(x)) - Fraction(frame['origin'][i])) / Fraction(frame['L'])
                       for i, x in enumerate(vertex)] for vertex in v], dtype=object)
        density /= frame['R']
    triangle = TriangleMesh(tuple(map(tuple, v.tolist())), tuple(map(tuple, f.tolist())))
    # Independent signed-tetrahedron implementation; one density-free integration
    # per region per frame. Density enters exactly once below.
    unit = compute_mass_properties(triangle, 1.)
    geometric = {'volume': unit.mass, 'first_moment': list(unit.first_moment),
                 'second_moment': [list(row) for row in unit.second_moment]}
    moments = [density * x for x in _ten(unit.mass, unit.first_moment, unit.second_moment)]
    properties = {'mass': moments[0], 'moments': moments,
                  'center_of_mass': list(unit.center_of_mass),
                  'central_inertia': [[density*x for x in row] for row in unit.central_inertia]}
    require(properties['mass'] > 0, 'density-weighted region mass underflow')
    canonical(properties)
    return geometric, properties


def _aggregate(records):
    require(bool(records), 'final material must have positive mass')
    moments = [math.fsum(p['moments'][i] for p in records) for i in range(10)]
    mass = scalar(moments[0], 'final mass', positive=True)
    anchor = max(records, key=lambda p: p['mass'])['center_of_mass']
    local = [math.fsum(p['mass'] * (p['center_of_mass'][i] - anchor[i])
                       for p in records) / mass for i in range(3)]
    center = [anchor[i] + local[i] for i in range(3)]
    inertia_terms = [[[] for _ in range(3)] for _ in range(3)]
    for p in records:
        delta = [p['center_of_mass'][i] - anchor[i] - local[i] for i in range(3)]
        squared = math.fsum(x*x for x in delta)
        for i in range(3):
            for j in range(3):
                inertia_terms[i][j].append(p['central_inertia'][i][j] + p['mass'] *
                    ((squared if i == j else 0.) - delta[i]*delta[j]))
    result = {'mass': mass, 'moments': moments, 'center_of_mass': center,
              'central_inertia': [[math.fsum(x) for x in row] for row in inertia_terms]}
    canonical(result)
    return result


def _claim(value):
    keys(value, ('frame', 'order', 'moments'), 'solver claim')
    require(value['frame'] == FRAME, 'solver claim must be in input-normalized frame')
    require(isinstance(value['order'], (list, tuple)) and tuple(value['order']) == ORDER,
            'solver claim has wrong moment order')
    require(isinstance(value['moments'], (tuple, list, np.ndarray)) and len(value['moments']) == 10,
            'solver claim must contain ten moments')
    moments = [scalar(x, 'solver moment') for x in value['moments']]
    require(moments[0] > 0, 'solver mass claim must be positive')
    return {'frame': FRAME, 'order': list(ORDER), 'moments': moments}


def _recheck_problem(problem, *, reuse_validation=False):
    require(type(reuse_validation) is bool, 'reuse_validation must be boolean')
    if reuse_validation:
        from .spin_problem import _reuse_constructed_problem_validation
        return _reuse_constructed_problem_validation(problem)
    fresh = SpinArtifactProblem(original=problem.original,
        fixed_regions=[(r.region_id, r.mesh, r.density) for r in problem.fixed_regions],
        design_regions=[(r.region_id, r.mesh, r.density) for r in problem.design_regions],
        length_unit=problem.length_unit, mass_unit=problem.mass_unit, target=problem.target,
        minimum_outer_wall=problem.minimum_outer_wall, validation_scope=problem.validation_scope)
    require(fresh.input_digest == problem.input_digest, 'locked problem identity mismatch')
    return fresh


def _check(problem, outer, results, material, carved, claim):
    require(outer == problem.original, 'outer canonical geometry bytes changed')
    require(not material.explicit_empty, 'final material must be nonempty')
    outer_region, material_region, carved_region = (m.exact_region() for m in (outer, material, carved))
    design_by_id = {r.region_id: r for r in problem.design_regions}
    decompositions, retained_native, carved_native, region_records = [], [], [], []
    world_records, normalized_records = [], []

    def integrate_record(region_id, kind, mesh, density):
        if mesh.explicit_empty:
            region_records.append({'id': region_id, 'kind': kind, 'density': density,
                                   'mesh': mesh.identity, 'explicit_empty': True,
                                   'geometric': None, 'physical': None,
                                   'normalized_geometric': None, 'normalized': None})
            return
        geometric, physical = _integrate(mesh, density)
        normalized_geo, normalized = _integrate(mesh, density, problem.normalization)
        world_records.append(physical)
        normalized_records.append(normalized)
        region_records.append({'id': region_id, 'kind': kind, 'density': density,
                               'mesh': mesh.identity, 'explicit_empty': False,
                               'geometric': geometric, 'physical': physical,
                               'normalized_geometric': normalized_geo, 'normalized': normalized})

    fixed_native = []
    for r in problem.fixed_regions:
        fixed_native.append(r.mesh.exact_region())
        integrate_record(r.region_id, 'fixed', r.mesh, r.density)
    boundary_facts = []
    for region_id, retained, removed in results:
        original_design = design_by_id[region_id]
        kept_region, removed_region = retained.exact_region(), removed.exact_region()
        decompositions.append({'id': region_id, 'partition': partition(
            original_design.mesh.exact_region(), [kept_region, removed_region])})
        retained_native.append(kept_region)
        carved_native.append(removed_region)
        boundary_facts.append({'id': region_id, 'retained': kept_region.facts, 'carved': removed_region.facts})
        integrate_record(region_id, 'retained', retained, original_design.density)
    material_partition = partition(material_region, fixed_native + retained_native)
    carved_partition = partition(carved_region, carved_native)
    if carved.explicit_empty:
        clearance = {'completed': True, 'applicable': False, 'at_least_minimum': True,
                     'state': 'NOT_APPLICABLE', 'squared_distance_exact': None}
    else:
        clearance = native().check_exact_boundary_clearance(outer_region, carved_region, problem.minimum_outer_wall)
        require(clearance['completed'], 'exact outer-wall distance could not be computed')
        require(clearance['applicable'] and clearance['at_least_minimum'], 'minimum outer wall violated')
    bounded_voids = material_region.facts['bounded_void_components']
    inter_cavity_web = ('NOT_APPLICABLE'
                        if type(bounded_voids) is int and bounded_voids < 2
                        else 'UNKNOWN')
    physical, normalized = _aggregate(world_records), _aggregate(normalized_records)
    differences = [abs(a-b) for a, b in zip(normalized['moments'], claim['moments'])]
    require(all(math.isfinite(x) and x <= TOLERANCE for x in differences),
            'independent normalized moment claim mismatch: ' + str(differences))
    # Unit-density combined geometry is a density-free report only. Actual mass
    # is the aggregate of the authenticated material-labelled regions above.
    combined_geo, _ = _integrate(material, 1.)
    wall_thickness_audit = {
        'outer_to_newly_carved': 'NOT_APPLICABLE' if carved.explicit_empty else 'VERIFIED',
        'same_component_feature': 'UNKNOWN',
        'inter_cavity_web': inter_cavity_web,
        'overall': 'UNKNOWN'}
    return {'schema': 'sealed-spin-artifact-v1', 'validation_scope': SCOPE,
            'geometry_status': 'MESH_MATERIAL_VALIDATED',
            'project_status': 'PARTIAL_BLOCKED', 'solver_target_feasibility': 'UNKNOWN',
            'build_source_attestation': 'NOT_VERIFIED', 'physical_validation': 'NOT_IN_SCOPE',
            'cad_validation': 'UNSUPPORTED_NOT_VERIFIED',
            'hydrostatic_validation': 'UNSUPPORTED_NOT_VERIFIED',
            'deformation_validation': 'UNSUPPORTED_NOT_VERIFIED',
            'original_evidence_admission': 'BLOCKED_UNKNOWN',
            'original_obj_source_identity': 'UNKNOWN_NOT_ATTESTED',
            'wall_thickness_audit': wall_thickness_audit,
            'same_component_feature_thickness': 'UNKNOWN',
            'inter_cavity_web_thickness': inter_cavity_web,
            'manufacturing_admission': 'BLOCKED_UNKNOWN_WALL_THICKNESS',
            'problem_digest': problem.input_digest, 'normalization': problem.normalization,
            'outer': outer.identity, 'final_material': material.identity, 'newly_carved': carved.identity,
            'design_results': [{'id': i, 'retained': r.identity, 'carved': c.identity} for i,r,c in results],
            'boundaries': {'original': outer_region.facts, 'material': material_region.facts,
                           'newly_carved': carved_region.facts, 'design_results': boundary_facts},
            'decompositions': decompositions, 'material_partition': material_partition,
            'carved_partition': carved_partition, 'clearance': clearance,
            'material_components': material_region.facts['material_components'],
            'bounded_void_components': material_region.facts['bounded_void_components'],
            'newly_carved_components': carved_region.facts['material_components'],
            'regions': region_records, 'combined_density_free_geometry': combined_geo,
            'physical': physical, 'normalized': normalized,
            'solver_claim': claim, 'moment_comparison': {'frame': FRAME, 'order': list(ORDER),
                'absolute_errors': differences, 'absolute_tolerance': TOLERANCE,
            'status': 'WITHIN_ENGINEERING_THRESHOLD'}}


_CREATED_ARTIFACTS = {}


def _artifact_state_digest(artifact):
    def mesh_state(mesh):
        if mesh.explicit_empty:
            return {'identity': mesh.identity, 'explicit_empty': True}
        region = mesh.exact_region()
        operands = getattr(region, 'construction_operands', ())
        return {'identity': mesh.identity, 'region_facts': region.facts,
                'construction_operands': [operand.facts for operand in operands]}
    problem = artifact.problem
    return digest({'problem_digest': problem.input_digest, 'problem_identity': problem.identity,
        'outer': mesh_state(artifact.outer),
        'design_results': [{'id': region_id, 'retained': mesh_state(retained),
                            'carved': mesh_state(carved)}
                           for region_id, retained, carved in artifact.design_results],
        'final_material': mesh_state(artifact.final_material),
        'newly_carved': mesh_state(artifact.newly_carved),
        'claim': artifact.solver_claim, 'facts': artifact.facts,
        'scientific_digest': artifact.scientific_digest,
        'evidence': artifact.evidence, 'evidence_digest': artifact.evidence_digest})


def _register_created_artifact(artifact):
    key = id(artifact)
    reference = weakref.ref(artifact, lambda ref, key=key:
        _CREATED_ARTIFACTS.pop(key, None) if _CREATED_ARTIFACTS.get(key, (None,))[0] is ref else None)
    _CREATED_ARTIFACTS[key] = (reference, _artifact_state_digest(artifact))


def _reuse_created_artifact_validation(artifact):
    """Reuse only this live artifact's completed _check result before save.

    Bundle readback still reconstructs and validates every payload independently.
    """
    entry = _CREATED_ARTIFACTS.get(id(artifact))
    if entry is None or entry[0]() is not artifact:
        return False
    try:
        return (implementation_observation() == artifact.evidence and
                _artifact_state_digest(artifact) == entry[1])
    except Exception:
        return False


@dataclass(frozen=True, init=False)
class SealedSpinArtifact:
    """Immutable recomputable observations; no caller-provided passed flags.

    Use create, then save_spin_artifact/load_spin_artifact. Reload requires the
    expected problem digest held by the caller's original run binding.
    """
    problem: SpinArtifactProblem
    outer: MeshSnapshot
    design_results: tuple
    final_material: MeshSnapshot
    newly_carved: MeshSnapshot
    _claim_json: str
    _facts_json: str
    _evidence_json: str
    scientific_digest: str
    evidence_digest: str

    def __init__(self):
        raise TypeError('use SealedSpinArtifact.create with actual geometry and solver claim')

    @classmethod
    def create(cls, problem, *, outer, design_results, final_material, newly_carved,
               solver_claim, validation_scope=SCOPE, _reuse_problem_validation=False):
        require(validation_scope == SCOPE, 'unsupported validation profile')
        require(isinstance(problem, SpinArtifactProblem), 'expected locked SpinArtifactProblem')
        before = implementation_observation()
        problem = _recheck_problem(problem, reuse_validation=_reuse_problem_validation)
        claim = _claim(solver_claim)
        require(isinstance(design_results, (tuple, list)), 'design results must be explicit list')
        results, seen = [], set()
        expected = {r.region_id for r in problem.design_regions}
        for row in design_results:
            require(isinstance(row, (tuple, list)) and len(row) == 3,
                    'design result must contain only id, retained mesh, carved mesh')
            region_id = row[0]
            require(isinstance(region_id, str) and region_id in expected and region_id not in seen,
                    'duplicate or unknown design ID')
            seen.add(region_id)
            results.append((region_id, MeshSnapshot(row[1]), MeshSnapshot(row[2])))
        require(seen == expected, 'missing design ID')
        outer, material, carved = (MeshSnapshot(m) for m in (outer, final_material, newly_carved))
        facts = _check(problem, outer, tuple(results), material, carved, claim)
        require(implementation_observation() == before, 'implementation file observation changed during checking')
        obj = object.__new__(cls)
        attrs = {'problem': problem, 'outer': outer, 'design_results': tuple(results),
                 'final_material': material, 'newly_carved': carved,
                 '_claim_json': canonical(claim), '_facts_json': canonical(facts),
                 '_evidence_json': canonical(before), 'scientific_digest': digest(facts),
                 'evidence_digest': digest({'scientific_digest': digest(facts), 'observation': before})}
        for name, value in attrs.items():
            object.__setattr__(obj, name, value)
        _register_created_artifact(obj)
        return obj

    @property
    def facts(self):
        return parse(self._facts_json)

    @property
    def evidence(self):
        return parse(self._evidence_json)

    @property
    def solver_claim(self):
        return parse(self._claim_json)

    def revalidate(self):
        require(implementation_observation() == self.evidence, 'stale implementation observation')
        fresh = type(self).create(self.problem, outer=self.outer, design_results=self.design_results,
            final_material=self.final_material, newly_carved=self.newly_carved, solver_claim=self.solver_claim)
        require(fresh.scientific_digest == self.scientific_digest and fresh.facts == self.facts,
                'recomputed scientific facts mismatch')
        require(fresh.evidence_digest == self.evidence_digest, 'recomputed evidence mismatch')
        return fresh
