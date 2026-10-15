"""Historical full-JSON versus streaming payload comparison, no acceptance grant.

Only this migration checker deliberately loads the old JSON into memory. Normal
archive readers must not do so. Original JSON and failure records are never edited.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
from make_it_stand.verification.source_trace_archive import read_archive


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def legacy_events(data):
    """Consume one in-memory historical dictionary, preserving every payload field."""
    passes=data['solve'].pop('passes')
    counts=[len(p['trace']) for p in passes]
    yield dict(kind='header',trace_counts=counts,payload=data)
    for index,p in enumerate(passes):
        trace=p.pop('trace')
        yield dict(kind='pass',pass_index=index,payload=p)
        for j,t in enumerate(trace):
            yield dict(kind='iteration',pass_index=index,trace_index=j,payload=t)


def compare(legacy,archives):
    if not archives:raise ValueError('NO_ARCHIVES')
    paths=[Path(legacy),*[Path(p) for p in archives]]
    before={str(p.resolve()):sha(p) for p in paths}
    data=json.loads(Path(legacy).read_text())
    missing=object();count=0
    readers=[read_archive(p) for p in archives]
    try:
        for rows in itertools.zip_longest(legacy_events(data),*readers,fillvalue=missing):
            if any(row is missing for row in rows):raise ValueError('RECORD_COUNT_MISMATCH')
            if any(row!=rows[0] for row in rows[1:]):
                raise ValueError('PAYLOAD_MISMATCH at record '+str(count))
            count+=1
    finally:
        for reader in readers:reader.close()
    for path,expected in before.items():
        if sha(path)!=expected:raise ValueError('FILES_CHANGED_DURING_COMPARE')
    return dict(scope='HISTORICAL_REPLAY_NOT_CURRENT_PAPER_ACCEPTANCE',equal_records=count,
        all_payloads_equal=True,files=before,
        archive_bytes_equal=len({before[str(Path(p).resolve())] for p in archives})==1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--legacy',type=Path,required=True)
    parser.add_argument('--archives',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=compare(args.legacy,args.archives)
    with args.output.open('x') as out:json.dump(result,out,indent=2,sort_keys=True,allow_nan=False)
    print('Compared complete payload records:',result['equal_records'],flush=True)
