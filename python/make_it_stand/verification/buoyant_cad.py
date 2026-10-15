"""Explicit saved A→B-to-CAD adapter; never optimizes or relabels a result.

The existing CAD exporter consumes the genuine nested plane seal. A separate
versioned attachment preserves the outer composition and source A identities.
Its CAD coordinates remain the existing normalized model coordinates; the
physical coordinate/density map is explicit, not silently assumed to be mm.
"""
from pathlib import Path
import subprocess

from .buoyant_material import audit_saved_buoyant_case
from .faster_material import _result, export_spin_it_faster_quadric_material_step
from .spin_bundle import load_spin_artifact
from .spin_common import canonical, digest, require


def export_saved_buoyant_step(case, *, bundle_root, run_index, stage,
                             occt_executable, output_directory):
    require(type(run_index) is int and run_index >= 0 and stage in ('A','B'),
            'invalid buoyant CAD run/stage selection')
    checked = audit_saved_buoyant_case(case, bundle_root=bundle_root)
    require(run_index < len(case['runs']) and checked['runs'][run_index]['qualified'],
            'buoyant CAD requires independently qualified saved material')
    selected = next((s for s in case['runs'][run_index]['stages'] if s['stage'] == stage), None)
    require(selected is not None and selected['status'] == 'CONVERGED',
            'requested buoyant CAD stage did not converge')
    outer_identity = selected['seal_identity']
    if stage == 'B':
        require(outer_identity['operation'] == 'FLOATING_ISLAND_BUOYANT_B_FROM_SEALED_A',
                'B CAD requires the original A→B composition identity')
        plane_identity = outer_identity['b_seal']
    else:
        plane_identity = outer_identity
    require(plane_identity['operation'] == 'SEAL_BUOYANT_PLANE_A',
            'saved buoyant stage lacks its genuine plane seal')
    root = Path(bundle_root).resolve(strict=True)
    artifact = load_spin_artifact(root/selected['relative_bundle'],
                                 expected_problem_digest=selected['problem_digest'])
    if stage == 'B':
        source = case['runs'][run_index]['stages'][0]
        a = load_spin_artifact(root/source['relative_bundle'],
                               expected_problem_digest=source['problem_digest'])
        require(outer_identity.get('status') == 'CONVERGED' and
                outer_identity.get('solver_status') == 'CONVERGED' and
                outer_identity.get('seal_error') is None and
                outer_identity.get('solver_result') == selected['solver_result'] and
                outer_identity.get('solver_result_digest') == digest(selected['solver_result']),
                'B CAD outer status/solver binding is inconsistent')
        expected = {'source_a_scientific_digest': a.scientific_digest,
                    'source_a_evidence_digest': a.evidence_digest,
                    'artifact_scientific_digest': artifact.scientific_digest,
                    'b_artifact_scientific_digest': artifact.scientific_digest,
                    'b_artifact_evidence_digest': artifact.evidence_digest}
        require(all(outer_identity.get(k) == v for k,v in expected.items()),
                'B CAD outer identity differs from loaded A/B artifacts')
        composition = outer_identity.get('composition_binding', {})
        require(composition.get('source_artifact_scientific_digest') == a.scientific_digest and
                composition.get('source_artifact_evidence_digest') == a.evidence_digest,
                'B CAD composition differs from loaded source A')
    result = _result('CONVERGED', plane_identity, artifact)
    # Rehydration must not change a single identity field or fabricate a seal.
    require(result.identity == plane_identity, 'saved plane seal is incomplete or inconsistent')
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=False)
    binding = {'schema':'saved-buoyant-cad-binding-v1', 'case_id':case['id'],
               'case_scientific_digest':case['scientific_hash'], 'run':run_index, 'stage':stage,
               'source_operation':outer_identity['operation'],
               'source_result_digest':digest(outer_identity),
               'source_a_scientific_digest':selected.get('source_a_scientific_digest'),
               'nested_plane_seal_digest':result.scientific_digest,
               'artifact_scientific_digest':artifact.scientific_digest,
               'normalization':artifact.problem.normalization,
               'length_unit':artifact.problem.length_unit,'mass_unit':artifact.problem.mass_unit,
               'cad_coordinates':'INPUT_NORMALIZED_MODEL',
               'scope':'SOLVER_MATERIAL_CAD_METHOD_NOT_MANUFACTURING_OR_AUTHOR_CASE'}
    (output/'request.json').write_text(canonical(binding)+'\n')
    try:
        cad = export_spin_it_faster_quadric_material_step(result,occt_executable=occt_executable,
            brep_path=output/'material.brep',step_path=output/'material.step')
    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
        (output/'failure.json').write_text(canonical({'type':type(error).__name__,'message':str(error)})+'\n')
        raise
    norm = artifact.problem.normalization
    report = {**binding, 'cad':cad, 'material_region_ids':cad['material_region_ids'],
              'physical_mass':cad['raw_mass_moments']['step'][0]*norm['R']*norm['L']**3,
              'maximum_normalized_material_moment_error':cad['maximum_normalized_material_moment_error']}
    (output/'buoyant-cad.json').write_text(canonical(report)+'\n')
    return report
