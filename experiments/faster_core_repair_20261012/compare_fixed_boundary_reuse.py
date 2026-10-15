"""Same native checker, same saved material: full versus certified fixed-face reuse."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin

HERE=Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    native=spin._native_module()
    config=json.loads((HERE/'core-benchmark-inputs-job-01/output/faces-12288.json').read_text())
    run=HERE/'layered-formal-benchmark-faces-12288-seed-0-job-01/output/run'
    request=json.loads((run/'configuration.json').read_text())['request']
    outer=np.asarray(config['request']['outer_vertices'])
    low=outer.min(axis=0); spans=outer.max(axis=0)-low
    center=low+spans*.5; length=float(max(spans))
    origin=[Fraction(float(x)) for x in center]
    target=[Fraction(float(x)) for x in request['target_origin']]
    axes=[[Fraction(float(x)) for x in row] for row in request['target_axes']]
    tokens=[]
    for point in config['request']['design_vertices']:
        model=[(Fraction(float(point[k]))-origin[k])/Fraction(length)-target[k] for k in range(3)]
        row=[sum(model[k]*axes[k][j] for k in range(3)) for j in range(3)]
        tokens.append([f'{x.numerator}/{x.denominator}' for x in row])
    started=time.monotonic()
    proof=native.verify_fixed_rational_boundary(tokens,np.asarray(config['request']['design_faces'],dtype=np.int32))
    preparation=time.monotonic()-started
    rows=[]
    ignored={'backend','algorithm','exact_pair_tests','reused_fixed_pair_tests','fixed_proof_identity'}
    for index in (1,73):
        path=run/f'evaluation-{index:06d}.json'; material=json.loads(path.read_text())
        vertices=material['variable_rational_vertices']; faces=np.asarray(material['variable_faces'],dtype=np.int32)
        samples={}; facts={}
        for mode in ('full','reuse'):
            timings=[]
            for sample in range(4):
                started=time.monotonic()
                value=(native.audit_exact_rational_triangle_boundary(vertices,faces,'OutwardPositive')
                    if mode=='full' else native.audit_rational_boundary_with_fixed_proof(vertices,faces,'OutwardPositive',proof))
                timings.append(time.monotonic()-started)
                if sample and value!=facts[mode]: raise ValueError('NONDETERMINISTIC_BOUNDARY_FACTS')
                facts[mode]=value
            samples[mode]=timings
        equal={k:v for k,v in facts['full'].items() if k not in ignored}=={
            k:v for k,v in facts['reuse'].items() if k not in ignored}
        row=dict(index=index,material_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            scientific_facts_equal=equal,samples_seconds=samples,
            full_predicate_pairs=facts['full']['exact_pair_tests'],
            computed_predicate_pairs=facts['reuse']['exact_pair_tests'],
            reused_predicate_pairs=facts['reuse']['reused_fixed_pair_tests'],
            both_valid=facts['full']['boundary_valid'] and facts['reuse']['boundary_valid'])
        rows.append(row)
    report=dict(rows=rows,passed=all(r['scientific_facts_equal'] and r['both_valid'] for r in rows),
        reference_preparation_seconds=preparation,reference_identity=proof.reference_identity,
        native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        boolean_calls=0,source_calls=0,core_admission=False,
        scope='Two saved actual 12288-face materials; one cold and three warm audits per mode')
    (output/'report.json').write_text(json.dumps(report,sort_keys=True,allow_nan=False))
    print(json.dumps(report),flush=True)
    if not report['passed']: raise SystemExit(1)


if __name__=='__main__': main(Path(sys.argv[1]))
