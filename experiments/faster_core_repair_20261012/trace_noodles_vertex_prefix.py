"""Trace the certified bad vertex in preserved prefixes; no Boolean or repair."""
import argparse
import hashlib
import json
from pathlib import Path
import gmpy2
import numpy as np
from noodles_vertex_link_witness import link_components


def main(prefix, witness_path, output):
    witness=json.loads(witness_path.read_text())['witnesses'][0]
    target=np.array([float(gmpy2.mpq(x)) for x in witness['exact_coordinate']])
    rows=[]
    for path in sorted(prefix.glob('domain-*.json')):
        packet=json.loads(path.read_text())
        vertices=np.asarray(packet['vertices'])
        matches=np.flatnonzero(np.all(vertices==target,axis=1)).tolist()
        stars=[dict(vertex=index,**link_components(packet['faces'],index)) for index in matches]
        rows.append(dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         exact_binary64_coordinate_matches=matches,stars=stars))
    rejected=[row for row in rows if any(not s['is_single_cycle'] for s in row['stars'])]
    report=dict(schema='preserved-prefix-vertex-link-trace-v1',target_binary64=target.tolist(),
        witness_sha256=hashlib.sha256(witness_path.read_bytes()).hexdigest(),rows=rows,
        first_recorded_bad_star=rejected[0]['file'] if rejected else None,
        limitation='Exact binary64 coordinate match only; no cross-state feature correspondence or entity repair.',
        geometry_admission=False,optimization_admission=False)
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(dict(first_recorded_bad_star=report['first_recorded_bad_star'],
        matched_prefixes=sum(bool(r['stars']) for r in rows))))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('prefix','witness','output'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();main(args.prefix,args.witness,args.output)
