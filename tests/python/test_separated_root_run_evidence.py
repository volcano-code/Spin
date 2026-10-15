"""Preserved numerical/material journal from the new explicit profile.

The first resource job failed at downstream record verification. These tests
only validate its preserved journal, not resource success or a saved bundle.
"""
import json
from pathlib import Path
import shutil

import pytest

from make_it_stand.verification.cherchi_run_evidence import verify_cherchi_run_evidence

RUN = Path(__file__).resolve().parents[2]/'experiments/faster_core_repair_20261012/separated-closeout-small-seed-0-job-01/output/run'


def test_new_explicit_profile_checks_preserved_same_call_run():
    result = verify_cherchi_run_evidence(RUN)
    assert result['record_integrity_valid']
    assert result['accepted_updates'] > 0


def test_new_profile_cannot_omit_or_change_fixed_proof(tmp_path):
    root = tmp_path/'run'
    shutil.copytree(RUN, root)
    config = json.loads((root/'configuration.json').read_text())
    config['fixed_boundary_proof_identity'] = '0'*64
    (root/'configuration.json').write_text(json.dumps(config))
    with pytest.raises(ValueError, match='FIXED_BOUNDARY_PROOF_IDENTITY_MISMATCH'):
        verify_cherchi_run_evidence(root)
