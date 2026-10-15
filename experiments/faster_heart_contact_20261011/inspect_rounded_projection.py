"""Read saved rational input and independently audit its binary64 derivative.

No coordinate repair, welding or solver admission occurs here. A valid rounded
mesh can still represent a different set; topology, equality and moment error
are separate observations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import gmpy2
import numpy as np

from rounded_input import exact_moments, publish


def inspect_region(native, rows, faces):
    from make_it_stand.verification.spin_problem import MeshSnapshot
    authority = MeshSnapshot.from_rationals(rows, faces)
    vertices, triangles = authority.arrays()
    exact = exact_moments(rows, faces)
    displayed = exact_moments(vertices, faces)
    errors = [abs(a-b) for a, b in zip(exact, displayed)]
    displacement = max(abs(gmpy2.mpq(float(v))-gmpy2.mpq(q))
        for row, original in zip(vertices, rows) for v, q in zip(row, original))
    audit = native.audit_exact_triangle_boundary(vertices, triangles, 'AnyNonzero')
    record = dict(authority_identity=authority.identity,
        binary64_boundary_audit=audit, maximum_coordinate_rounding_exact=str(displacement),
        moment_errors_exact=[str(x) for x in errors],
        moment_errors=[float(x) for x in errors],
        binary64_region_valid=False, represented_sets_equal=None, set_relation=None)
    try:
        floating = MeshSnapshot((vertices, triangles))
    except (ValueError, RuntimeError) as error:
        record['binary64_region_error'] = str(error)
        return record
    record['binary64_region_valid'] = True
    relation = native.relate_exact_material_regions(authority.exact_region(), floating.exact_region())
    record['set_relation'] = relation
    if relation['completed']:
        record['represented_sets_equal'] = relation['equal']
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--admission', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    import _make_it_stand as native
    before = time.monotonic()
    admission = json.loads((args.admission/'admission.json').read_text())
    if admission['status'] not in ('INPUT_ADMITTED_FEASIBILITY_UNKNOWN', 'INPUT_CAPACITY_EXCLUDED',
                                    'INPUT_GEOMETRY_BLOCKED'):
        raise ValueError('UNSUPPORTED_OR_INCOMPLETE_INPUT_ADMISSION')
    reports = {}
    for name in ('outer', 'design', 'fixed'):
        source = args.admission/(name+'.rational.json')
        expected = admission['files'][source.name]
        if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise ValueError('STALE_RATIONAL_INPUT: '+name)
        payload = json.loads(source.read_text())
        print(json.dumps(dict(stage=name, status='INSPECTING_BINARY64_PROJECTION')), flush=True)
        report = inspect_region(native, payload['vertices'], payload['faces'])
        report['authority_file_sha256'] = expected
        publish(args.output/(name+'.json'), report)
        reports[name] = report
    valid = all(r['binary64_region_valid'] for r in reports.values())
    exact = all(r['represented_sets_equal'] is True for r in reports.values())
    publish(args.output/'summary.json', dict(schema='rounded-heart-input-projection-audit-v1',
        problem_digest=admission['problem_digest'], all_binary64_regions_valid=valid,
        all_represented_sets_equal=exact, solver_admission_granted=False,
        simulation_admission_granted=False, optimization_invocations=0,
        elapsed_seconds=time.monotonic()-before,
        native_path=native.__file__, native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()))


if __name__ == '__main__':
    main()
