"""Contract for encoding the registered actual-state flower replay."""
from pathlib import Path
import subprocess
import hashlib
import json
import sys


def validate_replay_contract(row):
    expected = dict(material_profile='MATCHED_COMPARISON_V1', frames=288, fps=24,
        playback_rates=[1, .05], selected_phase=0, renders_actual_states=True,
        vibration_amplification=1, intervals=[[0, 8], [0, .2]])
    if any(row.get(key) != value for key, value in expected.items()):
        raise ValueError('MATCHED_REPLAY_CONTRACT_MISMATCH')
    return row['frames'], row['fps']


def encode_frames(frames, video):
    frames, video = Path(frames), Path(video)
    if video.exists():
        raise FileExistsError(video)
    expected = {f'frame-{i:04d}.png' for i in range(1, 289)}
    if {p.name for p in frames.glob('*.png')} != expected:
        raise ValueError('FRAME_SEQUENCE_INCOMPLETE_OR_EXTRA')
    labels = (
        "drawtext=text='Before / After - same material - actual MuJoCo states':x=20:y=18:fontsize=22:fontcolor=white:box=1:boxcolor=black@0.5,"
        "drawtext=text='1x playback':x=20:y=h-40:fontsize=22:fontcolor=white:enable='lt(t,8)',"
        "drawtext=text='0.05x playback - displacement not amplified':x=20:y=h-40:fontsize=22:fontcolor=white:enable='gte(t,8)'")
    subprocess.run(['ffmpeg','-v','error','-n','-threads','1','-framerate','24',
        '-start_number','1','-i',str(frames/'frame-%04d.png'),'-vf',labels,
        '-frames:v','288','-c:v','libx264','-threads','1','-pix_fmt','yuv420p',
        '-movflags','+faststart',str(video)],check=True)


def main(render, output):
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    record = json.loads((render/'render.json').read_text())
    validate_replay_contract(record)
    for name, digest in record['files'].items():
        if sha(render/name) != digest:
            raise ValueError('RENDER_FILE_CHANGED')
    output.mkdir(parents=True, exist_ok=False)
    video = output/'flower-matched-comparison.mp4'
    encode_frames(render/'frames', video)
    actual = json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-of','json',str(video)]))
    stream = actual['streams'][0]
    if stream['nb_frames'] != '288' or float(stream['duration']) != 12.:
        raise ValueError('ENCODED_VIDEO_DURATION_MISMATCH')
    (output/'video.json').write_text(json.dumps(dict(
        render_sha256=sha(render/'render.json'), video_sha256=sha(video),
        physics_sha256=record['physics_sha256'], material_sha256=record['material_sha256'],
        material_profile=record['material_profile'], frames=288, fps=24,
        playback_rates=[1,.05], selected_phase=0, scope='SAVED_MUJOCO_REPLAY_NOT_NEW_SIMULATION'),indent=2))


if __name__=='__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
