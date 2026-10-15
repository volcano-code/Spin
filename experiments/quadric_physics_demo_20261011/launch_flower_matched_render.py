"""Render existing admitted flower trajectories, without restarting any solver."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job


def main():
    base = HERE/'paper-flower-method-v2/mujoco-blender-v1'
    bundle = HERE/'paper-flower-method-v2/table-scene-v1'
    physics = base/'physics-job-01/output'
    blender = Path('/home/lin/.local/opt/blender-5.1.2-linux-x64/blender')
    files = dict(launcher=Path(__file__), renderer=HERE/'flower_mujoco_blender.py',
        material_policy=HERE/'flower_render_material_policy.py', blender=blender,
        summary=physics/'summary.json', material=bundle/'table-scene.json')
    for path in (bundle/'parts').glob('*.obj'):
        files['mesh_'+path.name] = path
    for path in physics.glob('*.npz'):
        files['trajectory_'+path.name] = path
    result = run_source_job(base/'render-matched-job-02',
        [str(blender), '-b', '-t', '4', '--python', str(files['renderer']), '--',
         '--bundle', str(bundle), '--physics', str(physics), '--output', '{job}/output',
         '--material-profile', 'MATCHED_COMPARISON_V1'], files)
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)
    rendered = base/'render-matched-job-02/output'
    encoder = HERE/'flower_matched_encode.py'
    encode_files = dict(launcher=Path(__file__), encoder=encoder,
                       record=rendered/'render.json')
    for path in rendered.rglob('*'):
        if path.is_file():
            encode_files[str(path)] = path
    encoded = run_source_job(base/'encode-matched-job-01',
        ['/home/lin/anaconda3/bin/python3','-B',str(encoder),str(rendered),'{job}/output'],
        encode_files)
    print(encoded.status, encoded.returncode, flush=True)
    if encoded.status != 'PROCESS_EXITED' or encoded.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
