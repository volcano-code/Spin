"""Finite classification permits guarded trials without rewriting old FD gates."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


def module():
    path = HERE/'build_heart_guarded_admission.py'
    assert path.is_file(), 'missing explicit branch-aware Heart admission'
    spec = importlib.util.spec_from_file_location('heart_guarded_admission', path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def evidence():
    load = lambda name: json.loads((HERE/name/'output/report.json').read_text())
    events = lambda name: [json.loads(line) for line in (HERE/name/'output/events.jsonl').read_text().splitlines()]
    return dict(preparation=load('heart-formal-mapping-prepare-job-01'),
        cap_numeric=load('heart-cherchi-heart-first-cap-failure-derivative-job-01'),
        cap_classification=load('heart-cherchi-heart-first-cap-failure-connection-job-02'),
        small_numeric=load('heart-cherchi-heart-first-small-face-rejection-derivative-job-01'),
        small_classification=load('heart-small-face-finite-classification-job-01'),
        reference=load('heart-projective-independent-reference-job-01'),
        cap_events=events('heart-cherchi-heart-first-cap-failure-derivative-job-01'),
        small_events=events('heart-cherchi-heart-first-small-face-rejection-derivative-job-01'))


def test_existing_evidence_supports_guarded_trials_not_neighborhood_certification():
    report = module().assess_finite_evidence(**evidence())
    assert report['finite_evidence_ready']
    assert report['old_five_point_gate_passed'] is False
    assert report['neighborhood_certification'] is False
    assert report['case_complete'] is False
    assert report['cap_classifications']['REAL_CONSTRUCTION_BRANCH_CHANGE'] == 12
    assert report['small_classifications']['REAL_CONSTRUCTION_BRANCH_CHANGE'] == 42


@pytest.mark.parametrize('fault', ['missed_branch', 'unsupported_direction', 'wrong_jacobian_not_rejected',
                                 'missing_sample', 'failed_preparation'])
def test_missing_or_inconsistent_heart_evidence_blocks_trials(fault):
    data = deepcopy(evidence())
    if fault == 'missed_branch':
        index = next(row['index'] for row in data['small_classification']['rows']
                     if row['classification'] == 'REAL_CONSTRUCTION_BRANCH_CHANGE')
        base = data['small_events'][0]['topology_hash']
        next(row for row in data['small_events'] if row['index'] == index)['topology_hash'] = base
    elif fault == 'unsupported_direction':
        data['reference']['direction'] = 1
    elif fault == 'wrong_jacobian_not_rejected':
        data['reference']['wrong_jacobian_status'] = 'SUPPORTED'
    elif fault == 'missing_sample':
        data['cap_classification']['rows'].pop()
    else:
        data['preparation']['passed'] = False
    with pytest.raises(ValueError):
        module().assess_finite_evidence(**data)
