"""Recover a numerical terminal without fabricating a completed run."""
import hashlib
import json
import shutil

import pytest

from make_it_stand.verification import cherchi_run_evidence as evidence
from test_cherchi_formal_run_evidence import RUN


def interrupted_prefix(tmp_path):
    root = tmp_path/'run'
    shutil.copytree(RUN, root)
    (root/'complete.json').unlink()
    (root/'result.json').unlink()
    (root/'terminal-material-reference.json').unlink()
    events = [json.loads(line) for line in (root/'events.jsonl').read_text().splitlines()]
    events = events[:next(i+1 for i, row in enumerate(events)
                         if row.get('kind') == 'NEWTON_EVENT' and
                         row['event'].get('type') == 'attempt_terminal' and
                         row['event'].get('reason') == 'converged')]
    (root/'events.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in events))
    return root, events


def test_saved_convergence_prefix_is_recoverable_but_not_complete(tmp_path):
    root, _ = interrupted_prefix(tmp_path)
    load = getattr(evidence, 'load_cherchi_terminal_prefix', None)
    assert callable(load), 'saved accepted terminal cannot yet be loaded for postprocessing'
    checkpoint = load(root)
    assert checkpoint['accepted_updates'] > 0
    assert checkpoint['terminal_event']['raw_equation_residual'] <= 1e-12
    assert checkpoint['terminal_event']['multiplier_mu'] > 0
    assert checkpoint['material']['quadric_coefficients'] == checkpoint['terminal_event']['quadric_coefficients']
    assert checkpoint['material_audit_complete'] is False
    assert checkpoint['run_complete'] is False
    with pytest.raises(ValueError):
        evidence.verify_cherchi_run_evidence(root)


@pytest.mark.parametrize('fault', ['tampered_material', 'missing_convergence', 'wrong_backend',
                                 'rehashed_stale_jacobian', 'no_updates'])
def test_invalid_terminal_prefix_is_rejected(tmp_path, fault):
    root, events = interrupted_prefix(tmp_path)
    load = getattr(evidence, 'load_cherchi_terminal_prefix', None)
    assert callable(load), 'terminal prefix verification is missing'
    saved = [row for row in events if row['kind'] == 'MAP_SAVED'][-1]
    if fault in ('tampered_material', 'rehashed_stale_jacobian'):
        path = root/saved['evaluation_file']
        material = json.loads(path.read_text())
        material['moment_jacobian'][0][0] += .001
        path.write_text(json.dumps(material))
        if fault == 'rehashed_stale_jacobian':
            saved['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    elif fault == 'missing_convergence':
        events.pop()
    elif fault == 'no_updates':
        for row in events:
            if row['kind'] == 'NEWTON_EVENT' and row['event']['type'] == 'evaluation':
                row['event']['step'] = 0.
    else:
        config = json.loads((root/'configuration.json').read_text())
        config['boolean_backend_identity'] = 'cgal'
        (root/'configuration.json').write_text(json.dumps(config))
    (root/'events.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in events))
    with pytest.raises(ValueError):
        load(root)
