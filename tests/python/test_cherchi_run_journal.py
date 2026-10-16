"""Evaluation checkpoints must make file entries durable before accepted events."""
import importlib
import importlib.util
import json
import os
from pathlib import Path
import pytest


def api():
    name='make_it_stand.verification.cherchi_run_journal'
    assert importlib.util.find_spec(name) is not None,'formal journal helper missing'
    return importlib.import_module(name)


def spy_fsync(monkeypatch):
    original=os.fsync;calls=[]
    def observe(fd):
        calls.append(os.readlink('/proc/self/fd/'+str(fd)))
        return original(fd)
    monkeypatch.setattr(os,'fsync',observe)
    return calls


def test_evaluation_checkpoint_preserves_records_and_orders_directory_sync(tmp_path,monkeypatch):
    root=tmp_path/'run';root.mkdir();calls=spy_fsync(monkeypatch)
    with api().RunJournal(root,durability_profile='EVALUATION_CHECKPOINT_V1') as journal:
        journal.publish('evaluation-000001.json',{'mass':1.})
        assert calls==[str(root/'evaluation-000001.json')]
        journal.emit({'kind':'STAGE_STARTED','stage':'integration'})
        journal.emit({'kind':'MAP_SAVED','evaluation':1})
        assert len(calls)==1
        assert len((root/'events.jsonl').read_text().splitlines())==2
        journal.emit({'kind':'NEWTON_EVENT','event':{'type':'evaluation','outcome':'accepted'}})
        assert calls[-2:]==[str(root),str(root/'events.jsonl')]
        assert journal.metrics['event_fsync_calls']==1
    assert json.loads((root/'evaluation-000001.json').read_text())=={'mass':1.}
    assert len((root/'events.jsonl').read_text().splitlines())==3


def test_default_still_syncs_each_event(tmp_path,monkeypatch):
    root=tmp_path/'run';root.mkdir();calls=spy_fsync(monkeypatch)
    with api().RunJournal(root) as journal:
        journal.emit({'kind':'STAGE_STARTED','stage':'source'})
        journal.emit({'kind':'STAGE_RETURNED','stage':'source'})
        assert calls==[str(root/'events.jsonl')]*2


def test_failure_close_syncs_pending_prefix(tmp_path,monkeypatch):
    root=tmp_path/'run';root.mkdir();calls=spy_fsync(monkeypatch)
    with pytest.raises(ValueError,match='registered failure'):
        with api().RunJournal(root,durability_profile='EVALUATION_CHECKPOINT_V1') as journal:
            journal.publish('evaluation-000001.json',{'mass':1.})
            journal.emit({'kind':'STAGE_FAILED','error':'registered failure'})
            raise ValueError('registered failure')
    assert calls[-2:]==[str(root),str(root/'events.jsonl')]
    assert json.loads((root/'events.jsonl').read_text())['kind']=='STAGE_FAILED'


def test_unknown_durability_is_rejected_before_opening_files(tmp_path):
    with pytest.raises(ValueError,match='unknown run durability profile'):
        api().RunJournal(tmp_path,durability_profile='UNREGISTERED')
    assert not (tmp_path/'events.jsonl').exists()

def test_failed_newton_event_preserves_infinite_diagnostic_without_aborting(tmp_path):
    with api().RunJournal(tmp_path) as journal:
        journal.emit({'kind':'NEWTON_EVENT','event':{'type':'evaluation','outcome':'rejected','raw_equation_residual':float('inf')}})
    row=json.loads((tmp_path/'events.jsonl').read_text())
    assert row['event']['raw_equation_residual']=={'$nonfinite_float':'+inf'}
    assert row['nonfinite_encoding']=='IEEE754_TAGGED_V1'

def test_nonfinite_accepted_material_still_rejected(tmp_path):
    with api().RunJournal(tmp_path) as journal:
        with pytest.raises(ValueError):
            journal.emit({'kind':'NEWTON_EVENT','event':{'type':'evaluation','outcome':'accepted','raw_equation_residual':float('inf')}})
        with pytest.raises(ValueError):journal.publish('evaluation-000001.json',{'mass':float('nan')})
