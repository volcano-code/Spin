"""Cross-case identity gates; these do not replace native artifact revalidation."""
from pathlib import Path
from teapot import STRATEGIES, check_files


def check_case_context(native, config, summary, strategy, problem_digest):
    from make_it_stand.verification.spin_common import digest
    if (strategy not in STRATEGIES or config['strategy']!=strategy or summary['strategy']!=strategy
        or native['solver_settings']['objective']!=STRATEGIES[strategy]
        or native['problem_digest']!=problem_digest or summary['problem_digest']!=problem_digest
        or digest(native)!=summary['scientific_digest'] or native['status']!=summary['status']
        or config['solver']!=native['solver_settings'] or config['geometry']!=native['source_settings']
        or config['boundary']!=native['boundary_settings'] or config['minimum_wall']!=native['minimum_outer_wall']
        or native['density']!=1. or config['meters_per_model_unit']!=.025 or config['density_kg_m3']!=1200.):
        raise ValueError('CASE_IDENTITY_OR_OBJECTIVE_MISMATCH')
    if summary['sealed'] and (native['status']!='BINARY_MOMENT_VALIDATED' or not native['artifact_scientific_digest']):
        raise ValueError('INVALID_SEALED_NUMERICAL_RESULT')


def check_collider_context(manifest, problem_digest):
    if manifest['input_digest']!=problem_digest:raise ValueError('COLLIDER_INPUT_MISMATCH')
    check_files(manifest['outer_hash']);check_files(manifest['source_hashes'])
    for p in manifest['profiles'].values():
        check_files(p['hashes'])
        if p['near_upright_passed'] is not True:raise ValueError('COLLIDER_SUPPORT_FAILED')


def check_video_context(report, problem_digest, trajectory):
    path=str(Path(trajectory).resolve())
    if (report['problem_digest']!=problem_digest or report['source_unchanged'] is not True
        or path not in report['trajectory_hashes']
        or not any(r.get('path')==path and r['status']=='RECORDED' for r in report['runs'])):
        raise ValueError('VIDEO_INPUT_OR_TRAJECTORY_MISMATCH')
    check_files(report['trajectory_hashes'])
