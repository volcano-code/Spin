import hashlib
import json

import pytest

from boundary_profile_inputs import select_saved_states


def fixture(tmp_path):
    rows = []
    for i in range(1, 5):
        name = f'evaluation-{i:06d}.json'
        payload = json.dumps({'quadric_coefficients': [i], 'geometry_hash': str(i)}).encode()
        (tmp_path/name).write_bytes(payload)
        rows.append(dict(kind='MAP_SAVED', evaluation=i, evaluation_file=name,
                         sha256=hashlib.sha256(payload).hexdigest()))
        rows.append(dict(kind='NEWTON_EVENT', event=dict(type='evaluation', evaluation_index=i,
            outcome='rejected' if i == 2 else 'accepted',
            reason='initial_evaluation_succeeded' if i == 1 else 'residual_decreased')))
    (tmp_path/'events.jsonl').write_text('\n'.join(json.dumps(x) for x in rows))
    (tmp_path/'terminal-material-reference.json').write_text(json.dumps(rows[6]))


def test_uses_first_actual_update_not_first_trial(tmp_path):
    fixture(tmp_path)
    rows = select_saved_states(tmp_path)
    assert [r['role'] for r in rows] == ['initial', 'first_accepted', 'terminal']
    assert [r['evaluation'] for r in rows] == [1, 3, 4]


def test_changed_saved_state_is_rejected(tmp_path):
    fixture(tmp_path)
    (tmp_path/'evaluation-000003.json').write_text('{}')
    with pytest.raises(ValueError, match='SAVED_STATE_CHANGED'):
        select_saved_states(tmp_path)
