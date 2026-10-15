"""Ordinary source build versus the frozen projective production materials."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value

HERE=Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    frozen=json.loads((HERE/'projective-formal-pipeline-job-01/output/registered-inputs.json').read_text())
    payload=dict(frozen['request'])
    payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    settings=spin.SpinEvaluatorConfig(**frozen['settings'])
    evaluator=InProcessCherchiEvaluator(request,settings,rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    native=spin._native_module()
    rows=[]
    for state in ('ellipsoid','one_sheet','two_sheet'):
        prior_path=HERE/'projective-formal-pipeline-job-01/output/fixed'/(state+'-material.json')
        prior=json.loads(prior_path.read_text())
        q=np.asarray(prior['quadric_coefficients'])
        start=dict(evaluator.call_counts)
        mapped=evaluator(q)
        current=evaluator.material(q)
        material_fields=('variable_vertices','variable_rational_vertices','variable_faces',
                         'quadric_source_vertices','quadric_source_faces','geometry_hash','topology_hash')
        numeric_fields=('variable_raw_moments','raw_moments','moment_jacobian')
        material_equal={field:current[field]==prior[field] for field in material_fields}
        numeric_equal={field:np.array_equal(mapped[field],prior[field]) for field in numeric_fields}
        cgal=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings)
        row=dict(state=state,prior_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest(),
            material_equal=material_equal,numeric_equal=numeric_equal,
            source_calls=evaluator.call_counts['source']-start['source'],
            boolean_calls=evaluator.call_counts['boolean']-start['boolean'],
            cgal_derivative_defined=cgal.derivative_defined,
            cgal_moment_error=float(np.max(np.abs(np.asarray(cgal.variable_raw_moments)-mapped['variable_raw_moments']))),
            cgal_jacobian_error=None if cgal.coefficient_jacobian is None else
                float(np.max(np.abs(np.asarray(cgal.coefficient_jacobian)-mapped['moment_jacobian']))))
        row['passed']=(all(material_equal.values()) and all(numeric_equal.values()) and
            row['source_calls']==row['boolean_calls']==1 and row['cgal_derivative_defined'] and
            row['cgal_moment_error']<=1e-10 and row['cgal_jacobian_error'] is not None and row['cgal_jacobian_error']<=1e-6)
        (output/(state+'-material.json')).write_text(json.dumps(current,default=_json_value,allow_nan=False))
        (output/(state+'-cgal.json')).write_text(json.dumps(asdict(cgal),default=_json_value,allow_nan=False))
        rows.append(row)
        (output/'rows.json').write_text(json.dumps(rows,allow_nan=False))
    report=dict(rows=rows,passed=all(r['passed'] for r in rows),core_admission=False,
        native=str(native.__file__),native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        scope='Three frozen projective production maps on ordinary source build; not a neighborhood or solver gate')
    (output/'report.json').write_text(json.dumps(report,allow_nan=False))
    if not report['passed']:
        raise ValueError('FORMAL_NATIVE_REBUILD_FIXED_MAP_MISMATCH')


if __name__=='__main__':main(Path(sys.argv[1]))
