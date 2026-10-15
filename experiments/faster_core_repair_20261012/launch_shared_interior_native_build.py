"""Independent ordinary build for exact shared-operand terminal construction."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    generation = sys.argv[1] if len(sys.argv) > 1 else '09'
    if generation not in ('09', '10', '11', '12'):
        raise ValueError('UNREGISTERED_SHARED_INTERIOR_BUILD')
    files = dict(launcher=Path(__file__), worker=HERE/'build_formal_native_core.py',
        cmake=ROOT/'CMakeLists.txt', boost_version=Path(sys.base_prefix)/'include/boost/version.hpp',
        regression=ROOT/'tests/python/test_exact_shared_interior_difference.py')
    if generation in ('11', '12'):
        files['classification_regression'] = ROOT/'tests/python/test_exact_material_boundary_classification.py'
    if generation == '12':
        files['replay_metadata_regression'] = ROOT/'tests/python/test_owned_region_replay_metadata.py'
    for folder in ('src', 'include', 'bindings', 'cmake'):
        for path in (ROOT/folder).rglob('*'):
            if path.is_file() and path.suffix in ('.cpp', '.h', '.hpp', '.cmake', '.txt'):
                files['source_'+str(path.relative_to(ROOT))] = path
    for name in ('libigl/CMakeLists.txt', 'libigl/cmake/LibiglDownloadExternal.cmake'):
        files[name] = ROOT/name
    result = run_serial_source_job(ROOT, HERE/('formal-native-core-build-job-'+generation),
        ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output', 'sharedvertex-slabs'], files)
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
