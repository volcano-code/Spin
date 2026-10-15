"""Authenticate durable preprocessing prefixes, never grant geometry admission."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_chain(stage_directories, fixed_parts, manifest_sha256, native_sha256,
                 *, require_complete=True):
    if not fixed_parts or len(set(fixed_parts)) != len(fixed_parts):
        raise ValueError('INVALID_FIXED_PART_LIST')
    excluded = []
    predecessor = None
    records = []
    for directory in stage_directories:
        directory = Path(directory)
        packet_path = directory / 'stage.json'
        packet = json.loads(packet_path.read_text())
        if packet.get('schema') != 'noodles-exact-domain-preprocessing-stage-v1':
            raise ValueError('STAGE_SCHEMA_MISMATCH')
        if packet.get('manifest_sha256') != manifest_sha256:
            raise ValueError('INPUT_IDENTITY_MISMATCH')
        if packet.get('native_sha256') != native_sha256:
            raise ValueError('NATIVE_IDENTITY_MISMATCH')
        if packet.get('predecessor_stage_sha256') != predecessor:
            raise ValueError('PREDECESSOR_CHAIN_MISMATCH')
        new = packet.get('excluded_parts', [])
        if (new != fixed_parts[:len(new)] or len(new) <= len(excluded)
                or len(new) > len(fixed_parts) or new[:len(excluded)] != excluded):
            raise ValueError('EXCLUSION_PREFIX_MISMATCH')
        if packet.get('part_count') != len(fixed_parts):
            raise ValueError('PART_COUNT_MISMATCH')
        complete = len(new) == len(fixed_parts)
        if packet.get('preprocessing_complete') is not complete:
            raise ValueError('FALSE_COMPLETION_STATUS')
        facts = packet.get('region_facts', {})
        if facts.get('completed') is not True or facts.get('valid') is not True:
            raise ValueError('INVALID_STAGE_REGION')
        if any(packet.get(key) is not False for key in
               ('geometry_admission', 'optimization_admission', 'simulation_admission')):
            raise ValueError('PREPROCESSING_CANNOT_GRANT_ADMISSION')
        boundary = directory / 'design.rational.json'
        if sha(boundary) != packet.get('design_sha256'):
            raise ValueError('BOUNDARY_IDENTITY_MISMATCH')
        predecessor = sha(packet_path)
        excluded = new
        records.append(dict(stage=str(packet_path.resolve()), sha256=predecessor,
                            boundary=str(boundary.resolve()), boundary_sha256=sha(boundary)))
    complete = len(excluded) == len(fixed_parts)
    if require_complete and not complete:
        raise ValueError('PREPROCESSING_INCOMPLETE')
    return dict(schema='noodles-preprocessing-chain-check-v1',
                preprocessing_complete=complete, excluded_part_count=len(excluded),
                expected_part_count=len(fixed_parts), stages=records,
                manifest_sha256=manifest_sha256, native_sha256=native_sha256,
                geometry_admission=False, optimization_admission=False,
                simulation_admission=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--stage', type=Path, action='append', required=True)
    parser.add_argument('--allow-prefix', action='store_true')
    args = parser.parse_args()
    manifest_path = args.input / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    # A chain must remain attached to all frozen source parts, not just its label.
    for name, entry in manifest['parts'].items():
        if sha(args.input / entry['model_obj']) != entry['sha256']['model']:
            raise ValueError('FROZEN_SOURCE_CHANGED:' + name)
    result = verify_chain(args.stage, manifest['construction']['fixed_parts'],
                          sha(manifest_path), sha(args.native),
                          require_complete=not args.allow_prefix)
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
