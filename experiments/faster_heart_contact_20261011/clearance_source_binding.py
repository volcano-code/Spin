"""Bind saved exact admission to the unchanged original Heart and fixed pin."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_source_contract(source, expected_spec):
    source=Path(source)
    if json.loads((source/'spec.json').read_text()) != expected_spec:
        raise ValueError('SOURCE_SPEC_REVISION_MISMATCH')
    manifest=json.loads((source/'preparation.json').read_text())['sources']
    for name,value in manifest.items():
        if sha(source/name) != value:raise ValueError('STALE_SOURCE: '+name)
    return {name:sha(source/name) for name in ('heart-body.obj','pin.obj','spec.json')}


def check_fixed_source_identity(source, candidate, authority_job_identity):
    for name,key in [('heart-body.obj','heart-body.obj'),('pin.obj','pin.obj'),('spec.json','spec')]:
        if candidate.get(name) != source[name]:
            raise ValueError('CANDIDATE_FIXED_SOURCE_MISMATCH: '+name)
        if authority_job_identity[key]['sha256'] != source[name]:
            raise ValueError('AUTHORITY_SOURCE_MISMATCH: '+name)


def check_candidate_authority_inputs(candidate, construction, authority_job_identity, expected_spec):
    """Reject altered fixed geometry before any expensive exact construction."""
    if construction['inherited_physics_and_solver'] != expected_spec:
        raise ValueError('INHERITED_SPEC_MISMATCH')
    frozen=construction['frozen_outer_source']
    for name,key in [('heart-body.obj','heart-body.obj'),('pin.obj','pin.obj'),('spec.json','spec')]:
        if frozen[name] != authority_job_identity[key]['sha256']:
            raise ValueError('AUTHORITY_SOURCE_MISMATCH: '+name)
        if name != 'spec.json' and sha(Path(candidate)/name) != frozen[name]:
            raise ValueError('CANDIDATE_FIXED_GEOMETRY_CHANGED: '+name)


def bind_saved_admission(source,candidate,old_authority,admission):
    from rounded_input import SPEC
    source,candidate,old_authority,admission=map(Path,(source,candidate,old_authority,admission))
    original=check_source_contract(source,SPEC)
    old_config=json.loads((old_authority.parent/'config.json').read_text())
    construction=json.loads((candidate/'construction.json').read_text())
    check_fixed_source_identity(original,construction['frozen_outer_source'],old_config['identity'])
    if construction['inherited_physics_and_solver'] != SPEC:
        raise ValueError('INHERITED_SPEC_MISMATCH')
    for name in ('heart-body.obj','pin.obj'):
        if sha(candidate/name) != original[name]:raise ValueError('CANDIDATE_FIXED_GEOMETRY_CHANGED')
    before=json.loads((old_authority/'admission.json').read_text())
    after=json.loads((admission/'admission.json').read_text())
    for folder,record in [(old_authority,before),(admission,after)]:
        for name,value in record['files'].items():
            if sha(folder/name) != value:raise ValueError('STALE_SAVED_ADMISSION: '+name)
    outer_sha=sha(old_authority/'outer.rational.json')
    if (after['outer_authority_source_sha256'] != outer_sha or
        sha(admission/'outer.rational.json') != outer_sha):
        raise ValueError('ORIGINAL_OUTER_AUTHORITY_CHANGED')
    if after['construction_sha256'] != sha(candidate/'construction.json'):
        raise ValueError('CONSTRUCTION_IDENTITY_MISMATCH')
    # Tie this geometric result to the payload actually executed in its envelope.
    new_config=json.loads((admission.parent/'config.json').read_text())
    for key,path in [('old-admission.json',old_authority/'admission.json'),
                     ('old-outer.rational.json',old_authority/'outer.rational.json'),
                     ('candidate-construction.json',candidate/'construction.json'),
                     ('candidate-pin.obj',candidate/'pin.obj'),
                     ('candidate-inner-source.obj',candidate/'inner-source.obj')]:
        if new_config['identity'][key]['sha256'] != sha(path):
            raise ValueError('EXECUTED_INPUT_IDENTITY_MISMATCH: '+key)
    return dict(schema='heart-clearance-source-binding-v1',source_identity_verified=True,
        frozen_source=original,outer_authority_sha256=outer_sha,
        admission_sha256=sha(admission/'admission.json'),
        construction_sha256=sha(candidate/'construction.json'),
        status=after['status'],optimization_invocations=0,simulation_invocations=0)


def main():
    parser=argparse.ArgumentParser()
    for name in ('source','candidate','old-authority','admission','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    result=bind_saved_admission(args.source,args.candidate,args.old_authority,args.admission)
    result['checker_sha256']=sha(__file__)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps(result))


if __name__ == '__main__':main()
