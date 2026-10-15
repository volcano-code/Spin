"""Force genuinely smaller steps: the first high-speed max-step cap was inactive."""
import json
from pathlib import Path
import argparse
import hashlib
import numpy as np
from heart_top import simulate
from run import write


def main():
    p=argparse.ArgumentParser(); p.add_argument('run',type=Path); a=p.parse_args()
    dest=a.run/'forced-step-refinement'; dest.mkdir(exist_ok=False)
    report=json.loads((a.run/'validation-complete/report.json').read_text())
    write(dest/'plan.json',dict(max_step_seconds=1/4000,hz=[35,45],direction=0,duration=8,
        reason='1/240 and1/480 caps inactive in high-speed adaptive integration; require actual evaluation-count change',
        angle_tolerance_degrees=.001))
    rows=[]
    for model in report['models']:
        for hz in [35,45]:
            name=model['model']; print(name,hz,flush=True)
            base=json.loads((a.run/f'supplemental-rotation/{name}-{hz}hz-dir0.json').read_text())
            fine=simulate(model['mass_kg'],model['center_m'],model['pivot_inertia_kg_m2'],
                          [0,0,2*np.pi*hz],max_step=1/4000)
            write(dest/f'{name}-{hz}hz-dir0.json',fine)
            error=float(np.max(np.abs(np.array(fine['tilt_degrees'])-base['tilt_degrees'])))
            changed=fine['function_evaluations']>1.5*base['function_evaluations']
            rows.append(dict(model=name,hz=hz,max_angle_difference_degrees=error,
                base_evaluations=base['function_evaluations'],fine_evaluations=fine['function_evaluations'],
                actual_integration_changed=changed,passed=changed and error<=.001 and
                max(fine['relative_energy_drift'],fine['relative_vertical_momentum_drift'])<=1e-6))
    write(dest/'report.json',dict(checks=rows,all_passed=all(r['passed'] for r in rows),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print(rows,flush=True)


if __name__=='__main__': main()
