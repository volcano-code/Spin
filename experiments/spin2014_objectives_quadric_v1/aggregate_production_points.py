"""Reintegrate 81 saved producer materials and apply the frozen FD protocol.

Run under source_job's bounded guardian, like point jobs. A terminal point
receipt alone never grants this gate. Large materials are streamed; only
coefficients, independently derived moments and branch metadata are retained.
"""
import argparse
from fractions import Fraction
from pathlib import Path
import time

from production_point_job import require_point_resource_envelope,save,start_disk_watch
from make_it_stand.verification.quadric_point_jobs import (
    _read,validate_point_manifest,verify_point_sources,iter_complete_point_jobs)
from make_it_stand.verification.surface_derivatives import compare_quadric_five_point_differential
from make_it_stand.verification.rational_text import fraction_from_rational_text as rational
from make_it_stand.verification.rational_text import fraction_to_rational_text as text


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--manifest-identity',required=True)
    parser.add_argument('--jobs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--job-envelope',type=Path)
    args=parser.parse_args()
    require_point_resource_envelope(args.job_envelope)
    manifest=_read(args.manifest);validate_point_manifest(manifest)
    if manifest['identity']!=args.manifest_identity:
        raise ValueError('wrong frozen production manifest')
    settings=manifest['mapping_settings'];paths=settings['source_paths']
    verify_point_sources(manifest,paths)
    if Path(paths['aggregator']).resolve()!=Path(__file__).resolve():
        raise ValueError('wrong loaded aggregator')
    from make_it_stand.verification import _spin2014_authority as checker
    from make_it_stand.verification import mass_moments,rational_text,quadric_point_jobs,surface_derivatives,quadric_branch
    from make_it_stand import spin
    native=spin._native_module()
    loaded={'native':native.__file__,'checker':checker.__file__,
            'mapping_interface':spin.__file__,'mass_moments':mass_moments.__file__,
            'rational_text':rational_text.__file__,'manifest_checker':quadric_point_jobs.__file__,
            'surface_derivatives':surface_derivatives.__file__,'quadric_branch':quadric_branch.__file__}
    if any(name not in paths or Path(path).resolve()!=Path(paths[name]).resolve() for name,path in loaded.items()):
        raise ValueError('wrong loaded aggregation dependencies')
    disk_done=start_disk_watch(args.output.parent)
    args.output.mkdir(exist_ok=False)
    save(args.output/'started.json',dict(manifest_sha256=manifest['identity'],
        status='STARTED_NOT_ACCEPTED',scope='PRODUCTION_DERIVATIVE_CHECK_NOT_OPTIMIZATION'))
    started=time.monotonic();records=[];base_jacobian=None
    branch_fields=('provenance_topology_hash','surface_type','effective_normal_subdivisions',
                   'source_truncation_radius','source_interface_face_count','source_cap_face_count',
                   'boolean_backend_identity')
    for packet in iter_complete_point_jobs(manifest,[args.jobs/f'{i:02d}' for i in range(81)]):
        index=packet['index'];value=packet['evaluation']
        record=dict(index=index)
        try:
            expected=native.normalize_spin2014_authority_coefficients(
                [text(Fraction(x)) for x in manifest['points'][index]['raw_coefficients']])
            if any(value[name]!=expected[name] for name in ('geometry_authority','dual_authority')):
                raise ValueError('saved coefficient authority differs from frozen actual parameters')
            independent=checker.check_independent_authority_material_mapping(value,
                expected_request_identity=settings['expected_request_identity'],
                expected_context_identity=settings['expected_context_identity'],
                arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            analytic=quadric_branch.exact_full_rank_invariants(value['geometry_authority'],
                arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            if analytic['subtype'] not in ('ELLIPSOID','ONE_SHEET_HYPERBOLOID','TWO_SHEET_HYPERBOLOID'):
                raise ValueError('point outside admitted regular full-rank analytic branch')
            record.update(geometry_authority=value['geometry_authority'],
                moments=independent['independent_total_moments_target'],
                branch={**{name:value[name] for name in branch_fields},
                        'analytic_subtype':analytic['subtype'],'analytic_signature':analytic['signature']},
                analytic_invariants=analytic,status='INDEPENDENTLY_REINTEGRATED')
            if index==0: base_jacobian=value['native_moment_jacobian']
        except Exception as error:
            record.update(status='FAILED',exception_type=type(error).__name__,exception=str(error))
        save(args.output/f'independent-{index:02d}.json',record)
        records.append(record)
        print(f'point {index}: {record["status"]}',flush=True)
        # The next iteration releases the previous large material packet.
    checks=[]
    for row in manifest['rows']:
        ids=row['indices'];column=row['column']
        checked=dict(row)
        if base_jacobian is None or any(records[i]['status']!='INDEPENDENTLY_REINTEGRATED' for i in ids):
            checked.update(passed=False,status='INDEPENDENT_FAILURE')
        else:
            result=compare_quadric_five_point_differential(
                [manifest['points'][i]['raw_coefficients'][column] for i in ids],manifest['base'][column],
                [records[i]['geometry_authority'] for i in ids],
                [records[i]['moments'] for i in ids],base_jacobian)
            result['weights']=[text(w) for w in result['weights']]
            fixed_branch=all(records[i]['branch']==records[0]['branch'] for i in ids)
            divisor=Fraction(manifest['points'][ids[-1]]['raw_coefficients'][column])-Fraction(
                manifest['points'][ids[0]]['raw_coefficients'][column])
            central_q=[(rational(a)-rational(b))/divisor for a,b in zip(
                records[ids[-1]]['geometry_authority'],records[ids[0]]['geometry_authority'])]
            central=[(rational(a)-rational(b))/divisor for a,b in zip(
                records[ids[-1]]['moments'],records[ids[0]]['moments'])]
            prediction=[sum(rational(base_jacobian[i][j])*central_q[j] for j in range(10)) for i in range(10)]
            checked.update(result,fixed_branch=fixed_branch,
                original_centered_per_moment_error_authority=[text(abs(a-b)/max(Fraction(1),abs(b)))
                                                              for a,b in zip(central,prediction)],
                passed=fixed_branch and result['within_tolerance'],
                status='FIXED_BRANCH' if fixed_branch else 'BRANCH_CHANGE_REQUIRES_EXACT_EVIDENCE')
        checks.append(checked)
    verify_point_sources(manifest,paths)
    if _read(args.manifest)!=manifest: raise ValueError('manifest changed during aggregation')
    passed=len(records)==81 and len(checks)==30 and all(row['passed'] for row in checks)
    save(args.output/'terminal.json',dict(schema='quadric-production-split-five-point-results-v1',
        manifest_sha256=manifest['identity'],production_points=81,rows=checks,
        independently_reintegrated=sum(r['status']=='INDEPENDENTLY_REINTEGRATED' for r in records),
        elapsed_seconds=time.monotonic()-started,full_derivative_gate=passed,
        original_failures_preserved=True,estimated_errors_subtracted=False,
        optimizer_invocations=0,goal_complete=False))
    disk_done.set()
    return 0 if passed else 1


if __name__=='__main__': raise SystemExit(main())
