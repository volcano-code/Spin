"""Re-read all registered actual points; local numerical and branch verdict."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))


def checked_rows(registration,records,jacobian):
    from make_it_stand.verification.surface_derivatives import compare_quadric_five_point_differential
    if len(records)!=81 or [r['index'] for r in records]!=list(range(81)):
        raise ValueError('ALL_81_UNIQUE_ORDERED_POINTS_REQUIRED')
    rows=[]
    for row in registration['rows']:
        ids=row['indices'];column=row['column']
        local=all(r['derivative_defined'] and max(r['moment_errors'])<=1e-6 and r['material_facts']['valid'] for r in (records[i] for i in ids))
        result=compare_quadric_five_point_differential(
            [registration['points'][i][column] for i in ids],registration['base'][column],
            [records[i]['normalized_coefficients'] for i in ids],
            [records[i]['independent_total_moments'] for i in ids],jacobian)
        fixed=all(records[i]['branch_identity']==records[0]['branch_identity'] for i in ids)
        rows.append(dict(column=column,factor=row['factor'],indices=ids,errors=result['per_moment_error'],
            numerical_pass=result['within_tolerance'],fixed_branch=fixed,
            passed=local and fixed and result['within_tolerance'],
            status='FIXED_BRANCH_CHECK' if fixed else 'BRANCH_CHANGE_REQUIRES_WITNESS'))
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    root=args.directory.resolve();path=root/'registration.json'
    registration=json.loads(path.read_text());digest=hashlib.sha256(path.read_bytes()).hexdigest()
    for file,expected in registration['sources'].items():
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=expected:raise ValueError('REGISTERED_SOURCE_CHANGED:'+file)
    records=[];jacobian=None
    for index in range(81):
        job=root/('point-job-'+str(index).zfill(2));receipt=json.loads((job/'result.json').read_text())
        if receipt['status']!='PROCESS_EXITED' or receipt['returncode']!=0:raise ValueError('INCOMPLETE_PRODUCTION_JOB:'+str(index))
        row=json.loads((job/'output/point.json').read_text())
        if row['registration_sha256']!=digest or row['index']!=index or row['raw_coefficients']!=registration['points'][index]:
            raise ValueError('PRODUCTION_POINT_IDENTITY_MISMATCH')
        evaluation=job/'output/actual-evaluation.json'
        if hashlib.sha256(evaluation.read_bytes()).hexdigest()!=row['evaluation_sha256']:raise ValueError('ACTUAL_EVALUATION_CHANGED')
        if index==0:jacobian=json.loads(evaluation.read_text())['coefficient_jacobian']
        records.append(row)
    if jacobian is None:raise ValueError('NO_PRODUCTION_BASE_JACOBIAN')
    rows=checked_rows(registration,records,jacobian)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(dict(schema='rounded-heart-faster-81-point-report-v1',
        registration_sha256=digest,production_points=81,rows=rows,
        derivative_admission=all(r['passed'] for r in rows),optimizer_invocations=0),indent=2)+'\n')


if __name__=='__main__':main()
