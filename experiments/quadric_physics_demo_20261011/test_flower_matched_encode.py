import importlib
import pytest
import json
import subprocess
from PIL import Image


def test_matched_encode_rejects_wrong_profile_or_playback():
    module = importlib.import_module('flower_matched_encode')
    row = dict(material_profile='MATCHED_COMPARISON_V1', frames=288, fps=24,
               playback_rates=[1, .05], selected_phase=0, renders_actual_states=True,
               vibration_amplification=1, intervals=[[0, 8], [0, .2]])
    assert module.validate_replay_contract(row) == (288, 24)
    for key, value in [('material_profile','PAPER_ROLES_V1'), ('frames',287),
                       ('playback_rates',[1,1]), ('selected_phase',3),
                       ('renders_actual_states',False), ('vibration_amplification',10),
                       ('intervals',[[0,8],[1,1.2]])]:
        with pytest.raises(ValueError):
            module.validate_replay_contract(dict(row, **{key:value}))


def test_encoder_produces_registered_duration_and_refuses_missing_frame(tmp_path):
    module = importlib.import_module('flower_matched_encode')
    frames = tmp_path/'frames'
    frames.mkdir()
    picture = Image.new('RGB', (32, 32), (80, 100, 120))
    for i in range(1, 289):
        picture.save(frames/f'frame-{i:04d}.png')
    video = tmp_path/'replay.mp4'
    module.encode_frames(frames, video)
    info = json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-of','json',str(video)]))
    assert info['streams'][0]['nb_frames'] == '288'
    assert float(info['streams'][0]['duration']) == 12.
    with pytest.raises(FileExistsError):
        module.encode_frames(frames, video)
    (frames/'frame-0042.png').unlink()
    with pytest.raises(ValueError, match='FRAME_SEQUENCE'):
        module.encode_frames(frames, tmp_path/'missing.mp4')
