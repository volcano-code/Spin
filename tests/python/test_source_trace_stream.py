"""Full diagnostic payload equivalence; no solver or physical success shortcut."""
import gc
import json
from dataclasses import asdict
import pytest
import _make_it_stand as native
from make_it_stand.verification import source_spin_material as bridge
from test_source_spin_material_bridge import problem, settings


def run(**changes):
    p=problem();cfg=settings(bridge)
    kw={key:asdict(value) for key,value in cfg.items()}
    kw['solver_settings'].update(changes)
    return native.solve_source_spin2014_material(*p.original.arrays(),
        *p.design_regions[0].mesh.arrays(),**kw,density=1.,minimum_outer_wall=.5)


def reconstruct(events):
    events=iter(events);header=next(events)
    assert header['kind']=='header'
    facts=header['payload'];facts['solve']['passes']=[]
    for event in events:
        if event['kind']=='pass':
            assert event['pass_index']==len(facts['solve']['passes'])
            p=event['payload'];p['trace']=[];facts['solve']['passes'].append(p)
        else:
            assert event['kind']=='iteration'
            p=facts['solve']['passes'][event['pass_index']]
            assert event['trace_index']==len(p['trace'])
            p['trace'].append(event['payload'])
    assert [len(p['trace']) for p in facts['solve']['passes']]==header['trace_counts']
    return facts


@pytest.mark.parametrize('changes',[{}, {'maximum_leaf_count':1},
    {'termination_policy':'PAPER_POST_TRANSITION_V1'}])
def test_stream_preserves_every_legacy_field_and_failure(changes,tmp_path):
    result=run(**changes)
    assert hasattr(result,'iter_diagnostics'), 'missing lazy full diagnostic iterator'
    expected=result.diagnostics
    # Serialization is allowed to change byte layout, never numeric facts.
    from make_it_stand.verification.source_trace_archive import write_archive,read_archive
    path=tmp_path/'native.jsonl'
    write_archive(path,result.iter_diagnostics())
    actual=reconstruct(read_archive(path))
    assert actual==expected


def test_iterator_owns_result_and_outputs_are_independent():
    result=run()
    assert hasattr(result,'iter_diagnostics'), 'missing lazy full diagnostic iterator'
    first=result.iter_diagnostics();second=result.iter_diagnostics()
    header=next(first);header['payload']['solve']['status']='FORGED'
    assert next(second)['payload']['solve']['status']=='BINARY_MOMENT_VALIDATED'
    del result,second;gc.collect()
    assert next(first)['kind']=='pass'
    remaining=list(first)
    assert sum(e['kind']=='iteration' for e in remaining)==39
    with pytest.raises(StopIteration):next(first)


@pytest.mark.parametrize('changes', [{}, {'maximum_leaf_count': 1}])
def test_encoded_native_events_preserve_all_payloads_and_archive_checks(changes, tmp_path):
    from make_it_stand.verification.source_trace_archive import write_encoded_archive, read_archive

    result = run(**changes)
    expected = list(result.iter_diagnostics())
    encoded = result.iter_diagnostics_json()
    del result
    gc.collect()
    path = tmp_path / 'encoded.jsonl'
    facts = write_encoded_archive(path, encoded)
    assert list(read_archive(path)) == expected
    assert facts['records'] == len(expected)
    with pytest.raises(FileExistsError):
        write_encoded_archive(path, [])


def test_trusted_native_writer_keeps_full_sequence_footer_and_science_hashes(tmp_path):
    from make_it_stand.verification.source_trace_archive import (
        read_archive, write_trusted_native_archive,
    )

    result = run()
    path = tmp_path / 'native-fast.jsonl'
    expected = list(result.iter_diagnostics())
    facts_header, facts = write_trusted_native_archive(path, result)
    actual = list(read_archive(path))

    assert actual == expected
    assert facts['records'] == 42
    assert facts_header == actual[0]
    assert facts_header['trace_counts'] == [38, 1]
    before = expected[0]['payload']['solve']
    after = actual[0]['payload']['solve']
    for key in ('input_hash', 'trace_hash', 'output_hash'):
        assert after[key] == before[key]
    footer = json.loads(path.read_bytes().splitlines()[-1])
    assert footer['kind'] == 'complete'
    assert footer['records'] == 42
    assert footer['payload_sha256'] == facts['payload_sha256']
    assert facts['archive_sha256'] == __import__('hashlib').sha256(path.read_bytes()).hexdigest()


def test_trusted_native_writer_rejects_caller_supplied_event_stream(tmp_path):
    from make_it_stand.verification.source_trace_archive import write_trusted_native_archive

    with pytest.raises(ValueError, match='NATIVE_RESULT_REQUIRED'):
        write_trusted_native_archive(tmp_path / 'forged.jsonl', object())
    assert not (tmp_path / 'forged.jsonl').exists()


def test_encoded_disk_reader_retains_signed_values_and_rejects_corruption(tmp_path):
    from test_source_online_trace import run as disk_run
    from make_it_stand.verification.source_trace_archive import write_encoded_archive, read_archive

    result = disk_run(tmp_path / 'trace')
    expected = [json.dumps(e, sort_keys=True, allow_nan=False) for e in result.iter_diagnostics()]
    path = tmp_path / 'disk.jsonl'
    write_encoded_archive(path, result.iter_diagnostics_json())
    assert [json.dumps(e, sort_keys=True, allow_nan=False) for e in read_archive(path)] == expected
    trace = next((tmp_path / 'trace').glob('*.slqptrace'))
    data = bytearray(trace.read_bytes()); data[100] ^= 0x40; trace.write_bytes(data)
    with pytest.raises(RuntimeError):
        write_encoded_archive(tmp_path / 'corrupt.jsonl', result.iter_diagnostics_json())
    with pytest.raises(ValueError):
        list(read_archive(tmp_path / 'corrupt.jsonl'))


def test_experiment_export_selects_explicit_encoding_and_retains_header(tmp_path):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'experiments/teapot_spin2014'))
    from online_compare import export_diagnostics
    from make_it_stand.verification.source_trace_archive import read_archive

    result = run()
    headers = []
    for encoding in ('dictionary', 'native-json'):
        path = tmp_path / (encoding + '.jsonl')
        header, facts = export_diagnostics(result, path, encoding)
        headers.append(header)
        events = list(read_archive(path))
        assert facts['records'] == len(events) == 42
        assert header == events[0]
    assert headers[0] == headers[1]
    with pytest.raises(ValueError, match='UNSUPPORTED_TRACE_ENCODING'):
        export_diagnostics(result, tmp_path / 'invalid.jsonl', 'skip-validation')
    assert not (tmp_path / 'invalid.jsonl').exists()
