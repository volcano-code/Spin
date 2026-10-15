"""Sequential, resource-bounded jobs for the approved reconstructed cases."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'python'))
NATIVE = ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
PYTHON = '/home/lin/anaconda3/bin/python3'


def launch_flower(folder, job, operation, application=None, *, material_directory=None,
                  saved_check=None, physics_directory=None, render_directory=None,delivery_directory=None):
    from make_it_stand.verification.source_job import run_source_job
    runner_root = ROOT/'experiments/spin_it_reproduction/faster'
    if operation == 'solve':
        worker = runner_root/'run_paper_table1.py'
        command = [str(worker), '--manifest', str(folder/'flower-cases.json'),
                   '--starts', '2', '--max-restarts', '8', '--output', '{job}/output/application.json']
        assets = [folder/'flower-cases.json', folder/'input.json', folder/'input-sha256.json']
    elif operation == 'check':
        worker = runner_root/'check_saved_planar_material.py'
        source=application or folder/'application.json'
        command = [str(worker), '--application', str(source), '--output', '{job}/output']
        assets = [source]
    elif operation == 'prepare':
        if application is None:
            raise ValueError('SAVED_APPLICATION_REQUIRED')
        worker = HERE/'paper_flower_output.py'
        command = [str(worker), '--folder', str(folder), '--application', str(application),
                   '--output', '{job}/output', '--shared-arrangement']
        assets = [application, *[p for p in folder.iterdir() if p.suffix in ('.json', '.obj')],
                  HERE/'paper_flower_fixture.py', HERE/'solid_view.py',
                  HERE/'paper_flower_shared_partition.py', HERE/'flower_grip_preview.py',
                  HERE/'flower-grip-preview-v1/assembly.json']
        assets += list((HERE/'flower-grip-preview-v1').glob('*.obj'))
    elif operation == 'qualify':
        worker = HERE/'paper_flower_output.py'
        command = [str(worker), '--exact-audit', '--output', str(folder),
                   '--audit-report-directory', '{job}/output']
        assets = list(folder.glob('*.obj'))+[folder/'output-preparation.json']
        assets += list(folder.glob('shared-arrangement.json'))
    elif operation == 'saved-check':
        if application is None:
            raise ValueError('SAVED_APPLICATION_REQUIRED')
        worker = HERE/'paper_flower_saved_check.py'
        derived = material_directory or folder/'prepare-job-01/output'
        digest = __import__('hashlib').sha256(NATIVE.read_bytes()).hexdigest()
        command = [str(worker), '--folder', str(folder), '--application', str(application),
            '--derived', str(derived), '--assembly-reference', str(HERE/'flower-grip-preview-v1/assembly.json'),
            '--expected-native-sha', digest, '--output', '{job}/output']
        assets = [application, *folder.glob('*.json'), *derived.glob('*.obj'),
            *derived.glob('*.json'), HERE/'flower-grip-preview-v1/assembly.json',
            HERE/'paper_flower_fixture.py', HERE/'paper_flower_shared_partition.py']
    elif operation in ('physics','render','display','deliver'):
        derived=material_directory or folder/'prepare-job-01/output'
        saved=saved_check or folder/'saved-independent-job-01/output/saved-independent-check.json'
        physics=physics_directory or folder/'physics-job-01/output'
        rendered=render_directory or folder/'render-job-01/output'
        if operation=='physics':
            worker=HERE/'flower_run_physics.py'
            command=[str(worker),'--folder',str(folder),'--material',str(derived),
                '--material-audit',str(saved.parent),'--saved-material-check',str(saved),
                '--grip',str(HERE/'flower-grip-preview-v1'),'--output','{job}/output']
        elif operation=='render':
            worker=HERE/'flower_render_physics.py'
            command=[str(worker),'--material',str(derived),'--physics',str(physics),
                '--paper-colors','--output','{job}/output']
        elif operation=='display':
            worker=HERE/'flower_paper_display.py'
            command=[str(worker),'--folder',str(folder),'--saved-check',str(saved),
                '--physics',str(physics),'--output','{job}/output',
                '--source-video',str(rendered/'flower-physical-comparison.mp4')]
            if delivery_directory is not None:command+=['--delivery-directory',str(delivery_directory)]
        else:
            worker=HERE/'flower_delivery.py'
            command=[str(worker),'--folder',str(folder),'--material',str(derived),
                '--saved-check',str(saved),'--output','{job}/output']
        assets=[*folder.glob('*.json'),folder/'solve-job-01/output/application.json',
            saved,*derived.glob('*.json'),*derived.glob('*.obj'),
            HERE/'flower-grip-preview-v1/assembly.json',HERE/'paper_flower_fixture.py',
            HERE/'flower_run_physics.py',HERE/'flower_grip_preview.py',
            HERE/'flower_render_physics.py',HERE/'flower_result_figure.py']
        if operation in ('render','display'):assets.extend(physics.glob('*'))
        if operation=='display':assets.extend(rendered.glob('*'))
        if operation=='display' and delivery_directory is not None:assets.extend(delivery_directory.glob('*.*'))
        for i in (1,2):
            assets.extend([folder/f'solve-job-0{i}/output/application.json',
                           folder/f'check-job-0{i}/output/scientific.json'])
    else:
        raise ValueError('UNREGISTERED_JOB_OPERATION')
    files = dict(controller=Path(__file__), worker=worker, loaded_native=NATIVE,
                 cpp_runtime=Path('/usr/lib/x86_64-linux-gnu/libstdc++.so.6'))
    for index, file in enumerate(assets):
        files['asset.'+str(index)+'.'+file.name] = file
    for name in ('spin.py', '__init__.py'):
        files['api.'+name] = ROOT/'python/make_it_stand'/name
    for name in ('planar_material.py', 'run_evidence.py', 'spin_common.py'):
        files['checker.'+name] = ROOT/'python/make_it_stand/verification'/name
    # The module is pinned by its actual bytes. This does not claim that a
    # current dirty source tree produced the module; A audits that separately.
    interpreter = '/home/lin/anaconda3/envs/geosplatting/bin/python' if operation == 'prepare' else PYTHON
    if operation == 'prepare':
        # Bind actual GEOS construction dependencies, not just their names.
        import subprocess
        discovery=subprocess.check_output([interpreter,'-c',
            'import json,shapely,shapely.lib,trimesh,numpy;print(json.dumps([shapely.__file__,shapely.lib.__file__,trimesh.__file__,numpy.__file__]))'],text=True)
        for index,path in enumerate(json.loads(discovery)):
            files['construction_dependency.'+str(index)]=Path(path)
        files['construction_shared_partition']=HERE/'paper_flower_shared_partition.py'
    result = run_source_job(job, ['/usr/bin/env',
        'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6', interpreter, '-B', *command], files,
        environment={'PYTHONPATH': str(ROOT/'python')+':'+str(NATIVE.parent)})
    receipt = dict(status=result.status, returncode=result.returncode, job=str(result.directory),
        elapsed_seconds=result.elapsed_ns/1e9, operation=operation,
        author_original=False, full_case_approval=False)
    print(json.dumps(receipt), flush=True)
    return result.returncode if result.returncode else (0 if result.status == 'PROCESS_EXITED' else 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('operation', choices=['solve', 'check', 'prepare', 'qualify', 'saved-check',
                                            'physics','render','display','deliver'])
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--job', type=Path, required=True)
    parser.add_argument('--application', type=Path)
    parser.add_argument('--material-directory', type=Path)
    parser.add_argument('--saved-check', type=Path)
    parser.add_argument('--physics-directory', type=Path)
    parser.add_argument('--render-directory', type=Path)
    parser.add_argument('--delivery-directory', type=Path)
    args = parser.parse_args()
    raise SystemExit(launch_flower(args.input.resolve(), args.job.resolve(), args.operation,
                                 args.application.resolve() if args.application else None,
                                 material_directory=args.material_directory,
                                 saved_check=args.saved_check,physics_directory=args.physics_directory,
                                 render_directory=args.render_directory,delivery_directory=args.delivery_directory))
