"""Offline symbolic reconstruction, without native objective/linear algebra calls.

Successful events currently omit the native bordered matrix: geodesic agreement
checks its action, not an elementwise comparison against a saved native matrix.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp


def formulas():
    r=sp.symbols('m lx ly lz qxx qyy qzz qyz qxz qxy')
    m,x,y,z,xx,yy,zz,yz,xz,xy=r
    a=yy+zz-(y*y+z*z)/m
    b=xx+zz-(x*x+z*z)/m
    c=xx+yy-(x*x+y*y)/m
    d=-xy+x*y/m
    f=(a*a+b*b+2*d*d)/(c*c)+z/m
    g=sp.Matrix([x,y,-xz+x*z/m,-yz+y*z/m])
    return tuple(sp.lambdify(r,v,'numpy',cse=True) for v in
                 (f,sp.Matrix([sp.diff(f,t) for t in r]),sp.hessian(f,r),g,
                  g.jacobian(r),[sp.hessian(v,r) for v in g]))


def dual(q):
    p=np.array(q,dtype=float,copy=True);p[7:]*=2
    return p/np.linalg.norm(p)


def main():
    here=Path(__file__).parent
    f,df,hf,g,dg,hg=formulas()
    reports=[]
    for seed in (0,1):
        path=here/f'small-top-inprocess-solver-sequence-01/seed-{seed}-job-01/output/events.jsonl'
        events=[json.loads(line) for line in path.read_text().splitlines()]
        rows=[];first_step=None
        for index,event in enumerate(events):
            e=event.get('evaluation_payload')
            if not e: continue
            r=np.asarray(e['solver_raw_moments']);p=dual(event['quadric_coefficients'])
            grad=np.asarray(df(*r)).ravel();G=np.asarray(dg(*r));constraints=np.asarray(g(*r)).ravel()
            multipliers=np.linalg.lstsq(np.column_stack((-G.T,p)),-grad,rcond=None)[0]
            lam,mu=multipliers[:-1],multipliers[-1]
            optimality=grad+mu*p-G.T@lam
            residual=np.linalg.norm(np.r_[optimality,constraints])
            rows.append(dict(evaluation_index=event['evaluation_index'],
                coordinate_error=float(np.max(np.abs(p-e['solver_coordinates']))),
                objective_error=abs(float(f(*r))-event['objective']),
                mu_error=abs(float(mu)-event['multiplier_mu']),
                raw_residual_error=abs(float(residual)-event['raw_equation_residual'])))
            # Reconstruct the first accepted direction away from roundoff stationarity.
            if first_step is None and event['outcome']=='accepted':
                successors=[v for v in events[index+1:] if v['type']=='evaluation' and
                            v['restart']==event['restart'] and v['outcome']=='accepted']
                if successors:
                    nxt=successors[0];J=np.asarray(e['solver_moment_jacobian'])
                    H=np.asarray(hf(*r))-np.einsum('i,ijk->jk',lam,np.asarray(hg(*r)))
                    A=np.block([[H@J+mu*np.eye(10),-G.T,p[:,None]],
                                [G@J,np.zeros((4,4)),np.zeros((4,1))],
                                [2*p[None,:],np.zeros((1,4)),np.zeros((1,1))]])
                    rhs=np.r_[-optimality,-constraints,0.]
                    solution=np.linalg.solve(A,rhs)
                    tangent=solution[:10]-p*np.dot(p,solution[:10]);step=nxt['step']
                    predicted=np.cos(step)*p+np.sin(step)*tangent/np.linalg.norm(tangent)
                    predicted/=np.linalg.norm(predicted)
                    first_step=dict(from_evaluation=event['evaluation_index'],to_evaluation=nxt['evaluation_index'],
                        matrix_condition=float(np.linalg.cond(A)),
                        linear_residual=float(np.max(np.abs(A@solution-rhs))),
                        geodesic_coordinate_error=float(np.max(np.abs(predicted-dual(nxt['quadric_coefficients'])))))
        maxima={k:max(row[k] for row in rows) for k in rows[0] if k!='evaluation_index'}
        passed=all(v<=1e-10 for v in maxima.values()) and first_step['geodesic_coordinate_error']<=1e-6
        reports.append(dict(seed=seed,evaluation_count=len(rows),maxima=maxima,first_step=first_step,
                            passed=bool(passed),evaluations=rows))
    report=dict(checks=reports,passed=all(r['passed'] for r in reports),
                scope='Independent symbolic objective, multipliers, residual and first Newton geodesic; native matrix capture still pending')
    with (here/'independent-eq14-pilot-01.json').open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='checks'}))
    print(json.dumps([{k:v for k,v in r.items() if k!='evaluations'} for r in reports]))
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
