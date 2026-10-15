"""Explicit post-hoc follow-up to the original 2/6/15 Hz protocol.

35/45 Hz were selected from a local linearization of actual sealed inertia,
NOT from hidden successful nonlinear runs. Original failures remain intact.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from heart_top import integrate,properties,simulate
from upright_linearization import upright_matrix
from run import write


def main():
    p=argparse.ArgumentParser(); p.add_argument('run',type=Path); a=p.parse_args()
    dest=a.run/'supplemental-rotation'; dest.mkdir(exist_ok=False)
    config=json.loads((a.run/'config.json').read_text())
    path=a.run/'validation/actual_meshes.npz'; arrays=np.load(path)
    write(dest/'plan.json',dict(hz=[35.,45.],directions=list(range(8)),seconds=8,
        tilt_degrees=2,max_step_seconds=1/240,refined_max_step_seconds=1/480,
        reason='POST_HOC_PROTOCOL_FROM_LOCAL_LINEARIZATION_NOT_ORIGINAL_PREDECLARED_MATRIX',
        geometry_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    records=[]; linear=[]
    for name,prefix in [('original','outer'),('carved','material')]:
        s=integrate(arrays[prefix+'_v'],arrays[prefix+'_f'])
        m,c,I=properties(s*config['physical_density_kg_per_m3']*
                         config['meters_per_model_unit']**np.array([3,4,4,4,5,5,5,5,5,5]))
        if np.linalg.norm(c[:2])>1e-12*np.linalg.norm(c): raise ValueError('local oracle requires axial CoM')
        for hz in [2.,6.,15.,30.,35.,40.,45.]:
            values=np.linalg.eigvals(upright_matrix(m,c[2],I,2*np.pi*hz))
            linear.append(dict(model=name,hz=hz,eigenvalues=[[v.real,v.imag] for v in values],
                max_real_eigenvalue=float(values.real.max()),scope='LOCAL_LINEAR_NOT_GLOBAL_CERTIFICATE'))
        for hz in [35.,45.]:
            for direction in range(8):
                print(f'{name} {hz:g} Hz direction{direction}',flush=True)
                q=simulate(m,c,I,[0,0,2*np.pi*hz],azimuth=direction*np.pi/4)
                if max(q['relative_energy_drift'],q['relative_vertical_momentum_drift'])>1e-6:
                    raise ValueError('invariant failure')
                write(dest/f'{name}-{hz:g}hz-dir{direction}.json',q)
                r=dict(model=name,hz=hz,direction=direction,**{k:v for k,v in q.items()
                    if k not in ('time','quaternions','tilt_degrees')})
                if direction==0:
                    fine=simulate(m,c,I,[0,0,2*np.pi*hz],max_step=1/480)
                    error=float(np.max(np.abs(np.array(q['tilt_degrees'])-fine['tilt_degrees'])))
                    r['step_refinement_max_angle_difference_degrees']=error
                    write(dest/f'{name}-{hz:g}hz-refined.json',fine)
                    if error>1e-3: raise ValueError('refinement failure')
                records.append(r)
    write(dest/'report.json',dict(scope='POST_HOC_IDEAL_FIXED_PIVOT_NOT_PAPER_ORIGINAL',
        scenarios=records,local_linearizations=linear,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),Path(__file__).with_name('heart_top.py'),Path(__file__).with_name('upright_linearization.py')]}))
    print('Completed all32 supplemental scenarios and4 refinements',flush=True)


if __name__=='__main__': main()
