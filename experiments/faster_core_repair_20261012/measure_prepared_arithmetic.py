"""Frozen-module arithmetic reuse experiment; not a full solver benchmark.

Both paths compute all moments and coefficients for the same registered first
128 faces. This isolates boundary ownership; it is not the production moving
face selection and cannot establish a full-evaluation speedup.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import time
import numpy as np
import _cherchi_rational_moments as arithmetic

HERE = Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    baseline_path = HERE / 'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    spec = importlib.util.spec_from_file_location('_baseline._cherchi_rational_moments', baseline_path)
    old_arithmetic = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old_arithmetic)
    rows = []
    for label in ('small', 'faces-12288'):
        source = HERE / f'layered-formal-benchmark-{label}-seed-0-job-06/output/run/evaluation-000001.json'
        raw = source.read_bytes()
        material = json.loads(raw)
        vertices = material['variable_rational_vertices']
        faces = np.asarray(material['variable_faces'], dtype=np.int64)
        selection = list(range(min(128, len(faces))))
        baseline = None
        times = {'legacy': [], 'prepared': []}
        # First observation per path and three same-process repeats. Alternate
        # order to disclose, rather than hide, order/cache sensitivity.
        for repeat in range(4):
            for mode in (('legacy', 'prepared') if repeat % 2 == 0 else ('prepared', 'legacy')):
                start = time.perf_counter_ns()
                if mode == 'legacy':
                    moments = old_arithmetic.integrate_moments(vertices, faces, 1.)
                    coefficients = old_arithmetic.transport_coefficients_selected(vertices, faces, 1., selection)
                else:
                    boundary = arithmetic.PreparedBoundary(vertices, faces)
                    moments = boundary.integrate_moments(1.)
                    coefficients = boundary.transport_coefficients_selected(1., selection)
                times[mode].append((time.perf_counter_ns() - start) / 1e9)
                values = (moments, coefficients)
                if baseline is None:
                    baseline = tuple(value.copy() for value in values)
                for actual, expected in zip(values, baseline):
                    np.testing.assert_array_equal(actual, expected)
        rows.append(dict(label=label, input_sha256=hashlib.sha256(raw).hexdigest(),
            faces=len(faces), selected_faces=selection, seconds=times,
            results_bitwise_equal=True))
        print(json.dumps(rows[-1]), flush=True)
    module = Path(arithmetic.__file__)
    (output / 'report.json').write_text(json.dumps(dict(
        rows=rows, module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
        baseline_module_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        scope='ARITHMETIC_REUSE_FIRST_128_FACES_NOT_PRODUCTION_SELECTION',
        full_evaluation_performance_passed=False, core_admission=False), indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
