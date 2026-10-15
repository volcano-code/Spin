"""Independent integration of saved material under target-frame unit change.

Not a new Boolean or Cherchi derivative certificate; those gates stay false.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from independent_tetra_moments import integrate
from case_coordinate_contract import normalize_target_moments, normalize_target_quadric

HERE=Path(__file__).resolve().parent

def main(output):
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for job in ('repaired-heart-production-moments-job-01',
                'repaired-heart-exact-production-moments-job-01'):
        path=HERE/job/'output/point-000-material.json'
        raw=path.read_bytes(); saved=json.loads(raw)
        evaluation=saved['evaluation']
        vertices=np.asarray(evaluation['material_vertices'])
        faces=evaluation['material_faces']
        length=2.57
        normalized=[[Fraction(float(x))/Fraction(length) for x in row] for row in vertices]
        model_moments=integrate(vertices,faces,1.)
        normalized_moments=integrate(normalized,faces,1.)
        scaling_error=float(np.max(np.abs(normalized_moments-normalize_target_moments(model_moments,length,1.))))
        saved_error=float(np.max(np.abs(model_moments-evaluation['variable_raw_moments'])))
        q=np.asarray(saved['coefficients']); p,jac=normalize_target_quadric(q,length)
        def polynomial(c,v):
            x,y,z=v.T
            return np.array([np.ones(len(v)),x,y,z,x*x,y*y,z*z,2*y*z,2*x*z,2*x*y]).T@c
        divisor=np.linalg.norm(q*np.array([1.,length,length,length,*([length**2]*6)]))
        polynomial_error=float(np.max(np.abs(polynomial(p,vertices/length)-polynomial(q,vertices)/divisor)))
        if max(scaling_error,saved_error)>1e-10 or polynomial_error>1e-12:
            raise ValueError('SAVED_MATERIAL_COORDINATE_EQUIVALENCE_FAILED')
        rows.append(dict(source=str(path),sha256=hashlib.sha256(raw).hexdigest(),
            vertices=len(vertices),faces=len(faces),length=length,
            saved_integration_error=saved_error,moment_scaling_error=scaling_error,
            polynomial_error=polynomial_error,normalized_coefficients=p.tolist(),
            coefficient_coordinate_jacobian=jac.tolist()))
    (output/'report.json').write_text(json.dumps(dict(passed=True,rows=rows,
        scope='Saved target-frame material unit conversion only',
        cherchi_production_derivative_admission=False,case_admission=False),sort_keys=True))

if __name__=='__main__':main(Path(sys.argv[1]))
