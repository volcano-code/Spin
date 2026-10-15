"""Launch one bounded saved-material proof; refuses an existing job root."""
import argparse
import json
from pathlib import Path
import gmpy2
import gmpy2.gmpy2
from make_it_stand.verification.source_job import run_source_job


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    scripts = Path(__file__).resolve().parent
    evidence = root / '.superpowers/sdd/2026-10-07-spin2014-objectives-quadric'
    files = dict(driver=Path(__file__), proof=scripts/'audit_split_coplanar_events.py',
        envelope_checker=scripts/'production_point_job.py', status=scripts/'production_status.py',
        manifest=evidence/'split-production-manifest-04/manifest.json',
        gmpy2_package=Path(gmpy2.__file__), gmpy2_extension=Path(gmpy2.gmpy2.__file__))
    helpers = ('audit_production_coplanar_patches.py', 'diagnose_independent_primal_secants.py',
        'diagnose_independent_transport.py', 'diagnose_independent_full_chain.py',
        'diagnose_saved_cap_implicit_precision.py', 'validate_production_five_point_map.py')
    for name in helpers:
        files['helper_'+name] = evidence/name
    for name in ('coplanar_patch.py', 'quadric_point_jobs.py', 'rational_text.py',
                 'surface_derivatives.py', 'source_job.py', 'mass_moments.py'):
        files['checker_'+name] = root/'python/make_it_stand/verification'/name
    files['mapping_interface'] = root/'python/make_it_stand/spin.py'
    manifest = json.loads(files['manifest'].read_text())
    for name, path in manifest['mapping_settings']['source_paths'].items():
        files['manifest_'+name] = Path(path)
    for index in [0, *range(57, 65), *range(73, 81)]:
        point = evidence/f'split-production-point-04-{index:03d}'
        for name in ('config.json', 'result.json', 'complete.json', 'kernel.json'):
            files[f'point_{index}_{name}'] = point/name
        for name in ('terminal.json', 'point.json', 'independent.json'):
            files[f'point_{index}_{name}'] = point/'output'/name
    result = run_source_job(args.job,
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3', '-B', str(files['proof']),
         '--evidence', str(evidence), '--output', '{job}/output',
         '--job-envelope', '{job}/kernel.json'], files,
        environment={'PYTHONPATH': str(root/'python')})
    print(json.dumps(dict(status=result.status, returncode=result.returncode,
        scope='EXACT_SAVED_EVENT_WITNESS_NOT_DERIVATIVE_ADMISSION')), flush=True)
    return 0 if result.returncode == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
