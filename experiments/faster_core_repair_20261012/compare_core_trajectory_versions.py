"""Cross-version numeric trajectory comparison, not same-version repetition."""
import hashlib
import json
from pathlib import Path
import re
import sys

from make_it_stand.verification.cherchi_run_evidence import scientific_projection


def version_projection(value):
    value=scientific_projection(value)
    def visit(node):
        if isinstance(node,dict):
            return {key:visit(item) for key,item in node.items() if key!='evaluation_map_identity'}
        if isinstance(node,list):return [visit(item) for item in node]
        if isinstance(node,str):return re.sub(r'evaluation_map=[0-9a-f]{64}',
            'evaluation_map=REGISTERED_BUILD_VARIANT',node)
        return node
    return visit(value)


def main(before,after,output):
    old=json.loads((before/'output/run/result.json').read_text())['identity']['native_run']
    new=json.loads((after/'output/run/result.json').read_text())['identity']['native_run']
    a=json.loads((before/'output/report.json').read_text())
    b=json.loads((after/'output/report.json').read_text())
    report=dict(numeric_events_equal=version_projection(old['events'])==version_projection(new['events']),
        complete_trace_equal=old['trace']==new['trace'],
        terminal_coefficients_equal=old['quadric_coefficients']==new['quadric_coefficients'],
        terminal_material_equal=a['material_digest']==b['material_digest'],
        raw_residual_equal=a['raw_equation_residual']==b['raw_equation_residual'],
        before_normal_seconds=a['normal_solver_seconds'],after_normal_seconds=b['normal_solver_seconds'],
        before_result_sha256=hashlib.sha256((before/'output/run/result.json').read_bytes()).hexdigest(),
        after_result_sha256=hashlib.sha256((after/'output/run/result.json').read_bytes()).hexdigest(),
        scope='Cross-version numeric/geometry trajectory; cost fields and registered map identity differ',
        same_version_repeat=False,core_admission=False)
    report['passed']=all(report[key] for key in ('numeric_events_equal','complete_trace_equal',
        'terminal_coefficients_equal','terminal_material_equal','raw_residual_equal'))
    with output.open('x') as stream:json.dump(report,stream,sort_keys=True)
    print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]))
