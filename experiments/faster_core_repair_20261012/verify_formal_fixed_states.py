"""Same-operand formal Cherchi/CGAL regression; not solve admission."""
import json
import sys
from pathlib import Path
from dataclasses import asdict
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value

HERE = Path(__file__).resolve().parent


def main(output, frozen=None):
    output.mkdir(parents=True, exist_ok=False)
    if frozen is None:
        frozen = json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())
    request_data = dict(frozen['request'])
    if isinstance(request_data.get('settings'), dict):
        request_data['settings'] = spin.SpinSolverConfig(**request_data['settings'])
    request = spin.SpinItFasterRequest(**request_data)
    settings = spin.SpinEvaluatorConfig(**frozen['settings'])
    adapter = InProcessCherchiEvaluator(request, settings)
    native = spin._native_module()
    rows = []
    for name, coefficients in frozen['states'].items():
        q = np.asarray(coefficients, dtype=float)
        q /= np.linalg.norm(q)
        row = dict(state=name, core_admission=False)
        try:
            mapped = adapter(q)
            material = adapter.material(q)
            source = adapter.last['source']
            cgal = spin.evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=settings)
            # Exact comparison reuses the very same physical Boolean operands.
            construction = native.libigl_exact_boolean_construction(
                request.design_vertices, request.design_faces,
                material['quadric_source_vertices'], material['quadric_source_faces'],
                adapter.last['operation'])
            if not construction['succeeded']:
                raise ValueError('CGAL_EXACT_CONSTRUCTION_FAILED')
            a = native.make_exact_material_region_from_rationals(construction['vertices'], construction['faces'])
            b = native.make_exact_material_region_from_rationals(
                material['variable_rational_vertices'], material['variable_faces'])
            row.update(cgal_derivative_defined=cgal.derivative_defined,
                moment_error=float(np.max(np.abs(np.asarray(cgal.variable_raw_moments)-mapped['variable_raw_moments']))),
                jacobian_error=None if cgal.coefficient_jacobian is None else float(np.max(np.abs(np.asarray(cgal.coefficient_jacobian)-mapped['moment_jacobian']))),
                exact_relation=dict(native.relate_exact_material_regions(a,b)),
                exact_cgal_valid=a.facts['valid'], exact_cherchi_valid=b.facts['valid'],
                source_vertices=source['vertices'], source_faces=source['faces'],
                cherchi_call_counts=dict(adapter.call_counts))
            (output/(name+'-material.json')).write_text(json.dumps(material,default=_json_value,allow_nan=False))
            (output/(name+'-cgal.json')).write_text(json.dumps(asdict(cgal),default=_json_value,allow_nan=False))
        except Exception as exc:
            row.update(failure_type=type(exc).__name__,failure=str(exc))
        rows.append(row)
        (output/(name+'-report.json')).write_text(json.dumps(row,default=_json_value,allow_nan=False))
    (output/'report.json').write_text(json.dumps(dict(rows=rows,core_admission=False),default=_json_value,allow_nan=False))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
