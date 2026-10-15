"""Archive sequence, integrity and bounded-reading contracts."""
import importlib.util
import json
from pathlib import Path
import pytest


def api():
    assert importlib.util.find_spec('make_it_stand.verification.source_trace_archive'), 'missing stream archive'
    from make_it_stand.verification import source_trace_archive
    return source_trace_archive


def events():
    yield {'kind':'header','trace_counts':[1,0],'payload':{'schema':'test','solve':{'status':'UNKNOWN'}}}
    yield {'kind':'pass','pass_index':0,'payload':{'index':0,'status':'MAXIMUM_ITERATIONS'}}
    yield {'kind':'iteration','pass_index':0,'trace_index':0,'payload':{'iteration':0,'state':[.5]}}
    yield {'kind':'pass','pass_index':1,'payload':{'index':1,'status':'UNKNOWN'}}


def test_roundtrip_retains_failures_and_every_record(tmp_path):
    m=api();path=tmp_path/'run.jsonl'
    facts=m.write_archive(path,events())
    assert list(m.read_archive(path))==list(events())
    assert facts['records']==4
    assert len(facts['payload_sha256'])==64
    with pytest.raises(FileExistsError):m.write_archive(path,events())


def test_encoded_writer_preserves_raw_bytes_without_reencoding(tmp_path):
    m = api()
    raw = [(json.dumps(e, indent=None) + '\n').encode() for e in events()]
    path = tmp_path / 'encoded.jsonl'
    facts = m.write_encoded_archive(path, iter(raw))
    assert path.read_bytes().startswith(b''.join(raw))
    assert list(m.read_archive(path)) == list(events())
    assert facts['records'] == 4


@pytest.mark.parametrize('mutation', ['duplicate_key', 'nonfinite', 'overflow', 'no_newline',
                                     'multiline', 'str', 'missing', 'limit', 'footer'])
def test_encoded_writer_cannot_bypass_integrity(mutation, tmp_path):
    m = api()
    raw = [(json.dumps(e) + '\n').encode() for e in events()]
    limit = m.DEFAULT_RECORD_LIMIT
    if mutation == 'duplicate_key':
        raw[2] = raw[2].replace(b'"trace_index": 0', b'"trace_index": 0,"trace_index": 0')
    elif mutation == 'nonfinite':
        raw[2] = raw[2].replace(b'0.5', b'NaN')
    elif mutation == 'overflow':
        raw[2] = raw[2].replace(b'0.5', b'1e9999')
    elif mutation == 'no_newline': raw[2] = raw[2].rstrip()
    elif mutation == 'multiline': raw[2] += b'{}\n'
    elif mutation == 'str': raw[2] = raw[2].decode()
    elif mutation == 'missing': raw.pop(2)
    elif mutation == 'limit': limit = 16
    else: raw.append(b'{"kind":"complete"}\n')
    path = tmp_path / 'bad.jsonl'
    with pytest.raises((ValueError, TypeError)):
        m.write_encoded_archive(path, raw, max_record_bytes=limit)
    with pytest.raises(ValueError):
        list(m.read_archive(path))


def test_incomplete_writer_cannot_create_complete_archive(tmp_path):
    m=api();path=tmp_path/'run.jsonl'
    with pytest.raises(ValueError):m.write_archive(path,iter(list(events())[:-1]))
    assert path.exists()  # preserve interrupted evidence, no success footer
    with pytest.raises(ValueError):list(m.read_archive(path))


@pytest.mark.parametrize('mutation',['missing','duplicate','reorder','wrong_index','nonfinite'])
def test_invalid_sequence_is_rejected(tmp_path,mutation):
    m=api();seq=list(events())
    if mutation=='missing':seq.pop(2)
    elif mutation=='duplicate':seq.insert(3,seq[2])
    elif mutation=='reorder':seq[1],seq[2]=seq[2],seq[1]
    elif mutation=='wrong_index':seq[2]['trace_index']=2
    else:seq[2]['payload']['state']=[float('nan')]
    with pytest.raises(ValueError):m.write_archive(tmp_path/'bad.jsonl',seq)


@pytest.mark.parametrize('mutation',['truncate','tamper','trailing','duplicate_key'])
def test_corrupt_archive_never_verifies(tmp_path,mutation):
    m=api();path=tmp_path/'run.jsonl';m.write_archive(path,events())
    data=path.read_bytes()
    if mutation=='truncate':data=data.rsplit(b'\n',2)[0]+b'\n'
    elif mutation=='tamper':data=data.replace(b'0.5',b'0.6')
    elif mutation=='trailing':data+=b'{}\n'
    else:data=data.replace(b'"trace_index":0',b'"trace_index":0,"trace_index":0')
    path.write_bytes(data)
    with pytest.raises(ValueError):list(m.read_archive(path))


def test_reader_has_explicit_per_record_limit(tmp_path):
    m=api();path=tmp_path/'run.jsonl';m.write_archive(path,events())
    with pytest.raises(ValueError,match='RECORD_LIMIT'):list(m.read_archive(path,max_record_bytes=32))


def test_producer_cannot_mutate_promised_counts_after_yield(tmp_path):
    m=api()
    def incomplete():
        header={'kind':'header','trace_counts':[1],'payload':{}}
        yield header
        header['trace_counts'].clear()
    with pytest.raises(ValueError,match='INCOMPLETE_ARCHIVE'):
        m.write_archive(tmp_path/'bad.jsonl',incomplete())


def test_reader_header_mutation_cannot_change_verification(tmp_path):
    m=api();path=tmp_path/'run.jsonl';m.write_archive(path,events())
    reader=m.read_archive(path)
    header=next(reader);header['trace_counts'].clear()
    assert len(list(reader))==3
