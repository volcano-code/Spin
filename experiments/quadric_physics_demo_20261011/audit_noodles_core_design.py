"""One independent represented-boundary check; no Nef/manufacturing sweep."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
import _make_it_stand as native

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'experiments/faster_core_repair_20261012'))
from independent_tetra_moments import integrate


def main(input_folder,output):
    output.mkdir(parents=True,exist_ok=False)
    packet_path=input_folder/'domain-081.json';packet=json.loads(packet_path.read_text())
    manifest_path=HERE/'noodles-dense-v5/manifest.json';manifest=json.loads(manifest_path.read_text())
    if packet['manifest_sha256']!=hashlib.sha256(manifest_path.read_bytes()).hexdigest() or \
        packet['excluded_parts']!=manifest['construction']['fixed_parts'] or len(packet['excluded_parts'])!=81:
        raise ValueError('INCOMPLETE_OR_CHANGED_FIXED_PARTS')
    with (output/'events.jsonl').open('x') as events:
        def event(status,**values):
            events.write(json.dumps(dict(status=status,**values))+'\n');events.flush();os.fsync(events.fileno())
        event('INDEPENDENT_DESIGN_BOUNDARY_AUDIT_STARTED',faces=len(packet['faces']))
        start=time.monotonic()
        audit=native.audit_exact_triangle_boundary(packet['vertices'],packet['faces'],'AnyNonzero')
        with (output/'boundary-audit.json').open('x') as stream:json.dump(dict(audit),stream)
        event('INDEPENDENT_DESIGN_BOUNDARY_AUDIT_TERMINAL',passed=audit['boundary_valid'],seconds=time.monotonic()-start)
        if audit['boundary_valid']:
            event('INDEPENDENT_DESIGN_INTEGRATION_STARTED')
            moments=integrate(packet['vertices'],packet['faces'],1.)
            value=np.asarray(native.integrate_raw_mass_moments(np.asarray(packet['vertices']),np.asarray(packet['faces']),1.))
            error=float(np.max(np.abs(moments-value)))
            report=dict(input_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
                design_boundary_valid=True,design_moments= moments.tolist(),
                native_moment_max_error=error,design_integration_passed=error<=1e-10,
                geometry_admission=False,optimization_admission=False,
                pending=['actual outer union','actual fixed material','whole input partition','capacity check'])
            event('INDEPENDENT_DESIGN_INTEGRATION_TERMINAL',maximum_error=error)
        else:
            report=dict(design_boundary_valid=False,geometry_admission=False,optimization_admission=False,
                        input_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
                        blocker='DERIVED_DESIGN_BOUNDARY_REJECTED')
        with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.input,a.output)
