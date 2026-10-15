import importlib.util
import json
import sys
from pathlib import Path
import pytest
from make_it_stand.verification.source_trace_archive import write_archive

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/teapot_spin2014'))


def api():
    assert importlib.util.find_spec('compare_trace_replay'), 'missing full legacy/stream replay comparator'
    import compare_trace_replay
    return compare_trace_replay


def fixture():
    return {'schema':'literal','solve':{'status':'UNKNOWN','passes':[
        {'index':0,'status':'MAXIMUM_ITERATIONS','trace':[
            {'iteration':0,'state':[.5],'objective_trials':[{'accepted':False}]}]},
        {'index':1,'status':'UNKNOWN','trace':[]}]}}


def test_every_payload_field_compared_and_differences_rejected(tmp_path):
    m=api();old=tmp_path/'legacy.json';old.write_text(json.dumps(fixture()))
    stream=tmp_path/'run.jsonl';seq=list(m.legacy_events(fixture()))
    assert [r['kind'] for r in seq]==['header','pass','iteration','pass']
    assert seq[2]['payload']['objective_trials']==[{'accepted':False}]
    write_archive(stream,iter(seq))
    result=m.compare(old,[stream]);assert result['equal_records']==4
    assert result['all_payloads_equal']
    seq[2]['payload']['objective_trials'][0]['accepted']=True
    bad=tmp_path/'bad.jsonl';write_archive(bad,seq)
    with pytest.raises(ValueError,match='PAYLOAD_MISMATCH'):m.compare(old,[bad])


def test_no_archives_cannot_claim_equality(tmp_path):
    m=api();old=tmp_path/'legacy.json';old.write_text(json.dumps(fixture()))
    with pytest.raises(ValueError):m.compare(old,[])
