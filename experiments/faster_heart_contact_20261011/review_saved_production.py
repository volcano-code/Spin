"""Recheck a frozen production map with exact retriangulation evidence.

This reuses the registered integrator and five-point protocol. It never reruns
an optimizer or edits the old failed result. Its verdict is local to the named
manifest, NOT a permit for another input or another native mapping.
"""
import argparse
from fractions import Fraction
import hashlib
from pathlib import Path
import sys
import time

from make_it_stand.verification.quadric_point_jobs import (
    _read, validate_point_manifest, verify_point_sources, iter_complete_point_jobs)
from make_it_stand.verification.surface_derivatives import compare_quadric_five_point_differential
from make_it_stand.verification.rational_text import fraction_to_rational_text
from retriangulation_witness import compare_retriangulation


ROOT = Path(__file__).resolve().parents[2]


def checked_rows(manifest, records, jacobian):
    validate_point_manifest(manifest)
    if len(records) != 81 or any(type(r.get('index')) is not int or r['index'] != i
                                for i,r in enumerate(records)):
        raise ValueError('all 81 ordered independent records required')
    rows = []
    for row in manifest['rows']:
        ids, column = row['indices'], row['column']
        checked = dict(row)
        if jacobian is None or any(records[i]['status'] != 'INDEPENDENTLY_REINTEGRATED' for i in ids):
            checked.update(passed=False, status='INDEPENDENT_FAILURE')
        else:
            result = compare_quadric_five_point_differential(
                [manifest['points'][i]['raw_coefficients'][column] for i in ids], manifest['base'][column],
                [records[i]['geometry_authority'] for i in ids],
                [records[i]['moments'] for i in ids], jacobian)
            result['weights'] = [fraction_to_rational_text(x) for x in result['weights']]
            fixed = all(records[i]['branch'] == records[0]['branch'] for i in ids)
            # Only triangulation provenance may change under the exact witness.
            def without_triangulation(branch):
                return {k:v for k,v in branch.items() if k != 'provenance_topology_hash'}
            witnessed = all(records[i]['retriangulation'].get('proved') is True and
                without_triangulation(records[i]['branch']) == without_triangulation(records[0]['branch'])
                for i in ids)
            supported = fixed or witnessed
            checked.update(result, fixed_branch=fixed, exact_retriangulation_supported=witnessed,
                passed=supported and result['within_tolerance'],
                status=('FIXED_BRANCH' if fixed else 'EXACT_RETRIANGULATION_WITNESS' if witnessed else
                        'BRANCH_CHANGE_REQUIRES_EXACT_EVIDENCE'))
        rows.append(checked)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--manifest-identity',required=True)
    parser.add_argument('--jobs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--job-envelope',type=Path,required=True)
    args = parser.parse_args()
    sys.path.insert(0,str(ROOT/'experiments/spin2014_objectives_quadric_v1'))
    from production_point_job import require_point_resource_envelope, save, start_disk_watch
    require_point_resource_envelope(args.job_envelope)
    manifest = _read(args.manifest)
    validate_point_manifest(manifest)
    if manifest['identity'] != args.manifest_identity:
        raise ValueError('wrong frozen manifest')
    settings = manifest['mapping_settings']; paths = settings['source_paths']
    verify_point_sources(manifest,paths)
    from make_it_stand.verification import _spin2014_authority as checker
    from make_it_stand.verification import mass_moments,rational_text,quadric_point_jobs,surface_derivatives,quadric_branch
    from make_it_stand import spin
    native = spin._native_module()
    loaded = dict(native=native.__file__,checker=checker.__file__,mapping_interface=spin.__file__,
        mass_moments=mass_moments.__file__,rational_text=rational_text.__file__,
        manifest_checker=quadric_point_jobs.__file__,surface_derivatives=surface_derivatives.__file__,
        quadric_branch=quadric_branch.__file__)
    if any(Path(path).resolve() != Path(paths[name]).resolve() for name,path in loaded.items()):
        raise ValueError('wrong loaded registered dependency')
    args.output.mkdir(exist_ok=False)
    started = time.monotonic(); disk_done = start_disk_watch(args.output)
    save(args.output/'started.json',dict(status='STARTED_NOT_ACCEPTED',manifest_sha256=manifest['identity'],
        protocol='FIVE_POINT_FIXED_PARENT_REGION_RETRIANGULATION_V1',
        scope='REGISTERED_CUBOID_PRODUCTION_MAPPING_ONLY_NOT_ROUNDED_HEART'))
    records = []; base = None; jacobian = None
    branch_fields = ('provenance_topology_hash','surface_type','effective_normal_subdivisions',
        'source_truncation_radius','source_interface_face_count','source_cap_face_count','boolean_backend_identity')
    for packet in iter_complete_point_jobs(manifest,[args.jobs/f'{i:02d}' for i in range(81)]):
        i = packet['index']; value = packet['evaluation']; record = dict(index=i)
        try:
            expected = native.normalize_spin2014_authority_coefficients(
                [fraction_to_rational_text(Fraction(x)) for x in manifest['points'][i]['raw_coefficients']])
            if any(value[name] != expected[name] for name in ('geometry_authority','dual_authority')):
                raise ValueError('coefficient authority differs from frozen positions')
            independent = checker.check_independent_authority_material_mapping(value,
                expected_request_identity=settings['expected_request_identity'],
                expected_context_identity=settings['expected_context_identity'],
                arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            analytic = quadric_branch.exact_full_rank_invariants(value['geometry_authority'],
                arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            if analytic['subtype'] not in ('ELLIPSOID','ONE_SHEET_HYPERBOLOID','TWO_SHEET_HYPERBOLOID'):
                raise ValueError('unadmitted regular analytic branch')
            if i == 0: base = value; jacobian = value['native_moment_jacobian']
            record.update(status='INDEPENDENTLY_REINTEGRATED',geometry_authority=value['geometry_authority'],
                moments=independent['independent_total_moments_target'],
                branch={**{name:value[name] for name in branch_fields},
                        'analytic_subtype':analytic['subtype'],'analytic_signature':analytic['signature']},
                retriangulation=compare_retriangulation(base,value))
        except Exception as exc:
            record.update(status='FAILED',exception_type=type(exc).__name__,exception=str(exc))
        save(args.output/f'independent-{i:02d}.json',record)
        records.append(record)
        print(f'point {i}: {record["status"]}',flush=True)
    rows = checked_rows(manifest,records,jacobian)
    verify_point_sources(manifest,paths)
    if _read(args.manifest) != manifest: raise ValueError('manifest changed during review')
    passed = len(rows) == 30 and all(row['passed'] for row in rows)
    save(args.output/'terminal.json',dict(schema='saved-production-exact-retriangulation-review-v1',
        protocol='FIVE_POINT_FIXED_PARENT_REGION_RETRIANGULATION_V1',
        manifest_sha256=manifest['identity'], sources=manifest['sources'], rows=rows,
        evaluator_settings=settings['evaluator'],
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        witness_sha256=hashlib.sha256(Path(__file__).with_name('retriangulation_witness.py').read_bytes()).hexdigest(),
        full_derivative_gate=passed, production_points=81,
        independently_reintegrated=sum(r['status']=='INDEPENDENTLY_REINTEGRATED' for r in records),
        scope='NAMED_REGISTERED_MAPPING_ONLY_NOT_GENERAL_SOLVER_OR_ROUNDED_HEART',
        actual_new_material_evaluations=0, optimizer_invocations=0,
        original_failures_preserved=True, estimated_errors_subtracted=False,
        rounded_heart_admission_granted=False, goal_complete=False,
        elapsed_seconds=time.monotonic()-started))
    disk_done.set()
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
