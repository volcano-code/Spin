"""Strict data and local file observations shared by the hollow artifact path."""
from __future__ import annotations

import hashlib
import importlib
import json
import math
from pathlib import Path

import numpy as np

SCOPE = 'HOLLOW_ONLY_MESH_MATERIAL_V1'
ORDER = ('M', 'LX', 'LY', 'LZ', 'QXX', 'QYY', 'QZZ', 'QYZ', 'QXZ', 'QXY')
FRAME = 'INPUT_NORMALIZED_V1'
TOLERANCE = 1e-6


def require(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, expected, label):
    require(isinstance(value, dict) and set(value) == set(expected),
            label + ': missing or unknown fields')


def scalar(value, label, *, positive=False):
    require(isinstance(value, (int, float, np.integer, np.floating))
            and not isinstance(value, (bool, np.bool_)), label + ': expected real number')
    require(not isinstance(value, np.floating) or value.dtype.itemsize <= 8,
            label + ': precision beyond binary64 is unsupported')
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(label + ': unrepresentable number') from exc
    require(math.isfinite(result) and (not positive or result > 0.),
            label + ': expected finite' + (' positive' if positive else '') + ' number')
    if isinstance(value, (int, np.integer)):
        require(int(result) == int(value), label + ': lossy integer conversion')
    return result


def vector(value, label):
    require(isinstance(value, (list, tuple, np.ndarray)) and len(value) == 3,
            label + ': expected three coordinates')
    return [scalar(x, label) for x in value]


def text_label(value, label):
    require(isinstance(value, str) and bool(value.strip()) and '\x00' not in value,
            label + ': expected nonempty text')
    return value


def _finite_tree(value):
    if value is None or type(value) in (str, bool, int):
        return
    if type(value) is float:
        require(math.isfinite(value), 'JSON contains a nonfinite or overflowing number')
    elif type(value) in (list, tuple):
        for item in value:
            _finite_tree(item)
    elif type(value) is dict:
        require(all(type(k) is str for k in value), 'JSON object keys must be strings')
        for item in value.values():
            _finite_tree(item)
    else:
        raise ValueError('unsupported strict JSON value')


def canonical(value):
    _finite_tree(value)
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                      allow_nan=False)


def parse(payload):
    def pairs(items):
        out = {}
        for key, value in items:
            require(key not in out, 'duplicate JSON key: ' + key)
            out[key] = value
        return out
    def bad_constant(value):
        raise ValueError('nonfinite JSON number: ' + value)
    value = json.loads(payload, object_pairs_hook=pairs, parse_constant=bad_constant)
    _finite_tree(value)
    return value


def byte_digest(payload):
    return hashlib.sha256(payload).hexdigest()


def digest(value):
    return byte_digest(canonical(value).encode('ascii'))


def native():
    try:
        module = importlib.import_module('_make_it_stand')
    except ImportError as exc:
        raise RuntimeError('native exact geometry support is unavailable') from exc
    required = ('make_exact_material_region', 'empty_exact_material_region',
                'check_exact_material_partition', 'check_exact_boundary_clearance')
    if not all(hasattr(module, name) for name in required):
        raise RuntimeError('native exact geometry support is unavailable')
    return module


def implementation_observation():
    # Disk observations only: this is neither build attestation nor a loaded-memory hash.
    from .run_evidence import _record
    files = {name: Path(__file__).with_name(name + '.py') for name in
             ('spin_common', 'spin_problem', 'spin_artifact', 'spin_bundle',
              'mass_moments', 'run_evidence')}
    files['native'] = Path(native().__file__).resolve()
    files['spin_modes'] = Path(__file__).parent.parent / 'spin.py'
    return {'scope': 'OBSERVED_FILES_ON_DISK_ONLY',
            'build_source_attestation': 'NOT_VERIFIED',
            'files': {name: _record(path, label=name) for name, path in files.items()}}
