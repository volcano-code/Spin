import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import pytest


def test_replay_binds_packet_backend_and_literal_operands(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location('saved_boolean', root/'experiments/faster_core_repair_20261012/replay_saved_boolean.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    binary = tmp_path/'bridge.so'; binary.write_bytes(b'bridge')
    calls = []
    def boolean(*args):
        calls.append(args)
        return dict(exact_vertices=[], faces=[])
    monkeypatch.setitem(sys.modules, '_cherchi_core_boolean', SimpleNamespace(__file__=str(binary), boolean=boolean))
    body = dict(context=dict(backend_identity='cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:'+hashlib.sha256(b'bridge').hexdigest()),
        design_vertices=[[1,2,3]], design_faces=[[0,0,0]], source_vertices=[[4,5,6]], source_faces=[[0,0,0]], operation='difference')
    body['operand_packet_sha256'] = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    packet = tmp_path/'packet.json'; packet.write_text(json.dumps(body))
    module.main(packet, tmp_path/'new-job/output.json')
    assert calls == [(body['design_vertices'], body['design_faces'], body['source_vertices'], body['source_faces'], 'difference')]
    body['source_vertices'][0][0] = 8
    packet.write_text(json.dumps(body))
    with pytest.raises(ValueError, match='CAPTURE_PACKET_CHANGED'):
        module.main(packet, tmp_path/'other.json')
    assert len(calls) == 1
