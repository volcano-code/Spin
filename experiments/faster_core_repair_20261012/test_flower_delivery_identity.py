"""Keep actual dynamics, complete scene coverage and full inertia bound."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
FLOWER = HERE.parent/'quadric_physics_demo_20261011/paper-flower-method-v2'


def api():
    path = HERE/'audit_flower_delivery_identity.py'
    assert path.is_file(), 'flower final delivery identity audit is missing'
    spec = importlib.util.spec_from_file_location('flower_identity', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def actual():
    return (json.loads((FLOWER/'mujoco-blender-v1/physics-job-01/output/summary.json').read_text()),
            json.loads((FLOWER/'table-scene-v1/physical-properties.json').read_text()))


def test_actual_flower_scene_matrix_uses_full_distinct_inertias():
    summary, properties = actual()
    report = api().check_physics_matrix(summary, properties)
    assert report['scene_count'] == 64
    assert report['full_inertia_bound']
    assert report['shared_support_parameters']


@pytest.mark.parametrize('fault', ['missing_scene', 'discard_cross_terms', 'default_inertia', 'separate_support'])
def test_incomplete_or_unfair_flower_physics_is_rejected(fault):
    summary, properties = map(deepcopy, actual())
    if fault == 'missing_scene':
        summary['scenes'].pop()
    elif fault == 'discard_cross_terms':
        summary['scenes'][0]['loaded_inertia'][0][1] = 0.
        summary['scenes'][0]['loaded_inertia'][1][0] = 0.
    elif fault == 'default_inertia':
        row = next(scene for scene in summary['scenes'] if scene['variant'] == 'optimized')
        row['loaded_inertia'] = summary['scenes'][0]['loaded_inertia']
    else:
        summary['scenes'][-1]['support']['frequency_hz'] = 12.
    with pytest.raises(ValueError):
        api().check_physics_matrix(summary, properties)
