"""No Boolean calls: classify the saved 81-point batch's edge changes."""
import argparse
import hashlib
import json
from pathlib import Path
from exact_planar_patches import patch_boundary_signature


def main(inputs,output):
    output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((inputs/'manifest.json').read_text())
    reference=None;rows=[];base_branch=None
    for index,q in enumerate(manifest['samples']):
        path=inputs/f'point-{index:03d}-material.json'
        data=json.loads(path.read_text())
        if data['actual_coefficients']!=q:raise ValueError('POINT_PARAMETER_MISMATCH')
        certificate=patch_boundary_signature(data['vertices'],data['faces'])
        if reference is None:reference=certificate['patches'];base_branch=data['source_branch']
        row=dict(index=index,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                 indexed_patch_boundaries_equal=certificate['patches']==reference,
                 source_branch_equal=data['source_branch']==base_branch,
                 certificate=certificate)
        rows.append(row)
        with (output/f'point-{index:03d}-certificate.json').open('x') as stream:json.dump(row,stream)
    complete=len(rows)==manifest['point_count']==81
    stable=all(row['indexed_patch_boundaries_equal'] for row in rows)
    branch=all(row['source_branch_equal'] for row in rows)
    report=dict(point_count=len(rows),complete=complete,
        exactly_coplanar_coalesced_indexed_boundaries_stable=stable,source_branch_stable=branch,
        local_retriangulation_classification=complete and stable and branch,
        full_branch_certification=False,
        remaining='Cross-state vertex provenance correspondence and embedding not certified by this checker',
        scope='Exact coplanar internal-edge witnesses only; never compares different-parameter entities for equality')
    with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();main(args.inputs,args.output)
