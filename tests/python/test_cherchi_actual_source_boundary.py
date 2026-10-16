"""Reject genuinely degenerate physical operands before entering Cherchi."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
import _make_it_stand as native
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator

def evaluator():
    e=object.__new__(InProcessCherchiEvaluator)
    e.native=native;e.design=np.zeros((4,3));e.request=SimpleNamespace(design_faces=[])
    e._stage=lambda name,callback:callback()
    e.bridge=SimpleNamespace(boolean=lambda *args: {'called':True})
    return e

def test_real_evaluation65_is_rejected_before_boolean():
    root=Path(__file__).resolve().parents[2]
    packet=json.loads((root/'experiments/faster_core_repair_20261012/inertia-boolean65-capture-job-07/output/capture/boolean-evaluation65.json').read_text())
    with pytest.raises(ValueError,match='ACTUAL_SOURCE_BOUNDARY_INVALID'):
        evaluator()._boolean_with_valid_source(packet['source_vertices'],packet['source_faces'],'intersection')

def test_legal_tiny_faces_are_not_rejected_by_area_threshold():
    vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]])*1e-12
    faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]]
    assert evaluator()._boolean_with_valid_source(vertices,faces,'intersection')['called']
