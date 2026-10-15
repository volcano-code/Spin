"""Frozen 81-point scheduling/accounting, never a scientific pass certificate.

Independent jobs may finish in any order. Their stored packets must all be
authenticated before the existing material/derivative checkers inspect them.
Nothing here substitutes metadata equality for independent geometry checks.
"""
import hashlib
import json
from pathlib import Path
import re

from .surface_derivatives import quadric_five_point_layout


def file_sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def _identity(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),
                                    allow_nan=False).encode()).hexdigest()


def freeze_point_manifest(base,steps, *, sources,mapping_settings):
    if not isinstance(sources,dict) or not {'input','native','runner','checker'} <= sources.keys() or \
       any(not isinstance(k,str) or not isinstance(v,str) or
           re.fullmatch('[0-9a-f]{64}',v) is None for k,v in sources.items()):
        raise ValueError('frozen source/build/checker identities required')
    if not isinstance(mapping_settings,dict) or \
       mapping_settings.get('coefficient_policy')!='COEFFICIENT_MPFR1024_V1' or \
       mapping_settings.get('precision_bits')!=1024 or mapping_settings.get('normal_level')!=2:
        raise ValueError('explicit frozen authority mapping settings required')
    trials,rows=quadric_five_point_layout(base,steps)
    value=dict(schema='quadric-production-point-manifest-v1',
               sources=dict(sources),mapping_settings=mapping_settings,
               base=list(trials[0]),steps=list(steps),
               points=[dict(index=i,raw_coefficients=list(raw)) for i,raw in enumerate(trials)],
               rows=[dict(column=c,factor=f,indices=list(ids)) for c,f,ids in rows],
               scientific_tolerance=1e-6,source_interval_width='1/'+str(2**384),
               production_queries=81,optimizer_invocations=0,scientific_gate_granted=False)
    # Round-trip detaches caller-owned settings and requires JSON-safe values.
    value=json.loads(json.dumps(value,allow_nan=False))
    return {**value,'identity':_identity(value)}


def validate_point_manifest(manifest):
    try:
        expected=freeze_point_manifest(manifest['base'],manifest['steps'],
                                       sources=manifest['sources'],
                                       mapping_settings=manifest['mapping_settings'])
    except (KeyError,TypeError) as error:
        raise ValueError('incomplete frozen point manifest') from error
    if manifest!=expected:
        raise ValueError('frozen point manifest changed or layout incomplete')


def verify_point_sources(manifest,paths):
    """Before/after job guard; hash actual files, not caller's pass flags."""
    validate_point_manifest(manifest)
    if set(paths)!=set(manifest['sources']):
        raise ValueError('incomplete production dependency set')
    for name,path in paths.items():
        if not Path(path).is_file() or file_sha256(path)!=manifest['sources'][name]:
            raise ValueError('STALE_RESULT: production dependency changed: '+name)


def _read(path):
    def pairs(items):
        value={}
        for key,item in items:
            if key in value: raise ValueError('duplicate saved JSON key')
            value[key]=item
        return value
    def constant(value):
        raise ValueError('nonfinite saved JSON value: '+value)
    try:
        with Path(path).open() as stream:
            return json.load(stream,object_pairs_hook=pairs,parse_constant=constant)
    except (OSError,json.JSONDecodeError) as error:
        raise ValueError('missing or truncated point job evidence') from error


def load_complete_point_jobs(manifest,directories):
    return list(iter_complete_point_jobs(manifest,directories))


def iter_complete_point_jobs(manifest,directories):
    """Authenticate saved actual packets only; callers still verify all science.

    Resource-cut prefixes and failed jobs cannot fill a complete 81-point set.
    No caller-supplied passed=True value is consumed or granted here.
    """
    validate_point_manifest(manifest)
    directories=list(directories)
    if len(directories)!=81: raise ValueError('all 81 unique completed jobs required')
    entries={}
    for directory in directories:
        directory=Path(directory)
        receipt=_read(directory/'terminal.json')
        if not isinstance(receipt,dict) or receipt.get('schema')!='quadric-point-job-terminal-v1' or \
           receipt.get('manifest_sha256')!=manifest['identity'] or receipt.get('status')!='COMPLETE' or \
           type(receipt.get('actual_material_evaluations')) is not int or \
           receipt['actual_material_evaluations']!=1:
            raise ValueError('incomplete, failed or stale point job')
        index=receipt.get('index')
        if type(index) is not int or not 0<=index<81 or index in entries:
            raise ValueError('duplicate or invalid point job identity')
        entries[index]=(directory,receipt)
    # Only small receipts are retained; material packets are opened one at a time.
    for index in range(81):
        directory,receipt=entries[index]
        path=directory/'point.json'
        if not path.is_file() or file_sha256(path)!=receipt.get('packet_sha256'):
            raise ValueError('saved production packet changed or missing')
        independent=directory/'independent.json'
        if not independent.is_file() or file_sha256(independent)!=receipt.get('independent_sha256'):
            raise ValueError('saved independent check changed or missing')
        packet=_read(path)
        if not isinstance(packet,dict) or type(packet.get('index')) is not int or \
           packet['index']!=index or packet.get('raw_coefficients')!=manifest['points'][index]['raw_coefficients'] or \
           'exception' in packet or not isinstance(packet.get('evaluation'),dict):
            raise ValueError('point packet failed or differs from frozen input')
        yield packet
