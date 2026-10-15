"""Account for the entire registered matrix; never infer physical success.

This consumes saved trajectories/failures, not seals or solver pass flags.
Time-step and proxy differences are measurements, not a universal tilt cutoff.
"""
import math
from contact_matrix import scenarios


METRICS = ('maximum_tilt_degrees', 'final_tilt_degrees', 'final_spin_hz', 'maximum_penetration_m')


def summarize(protocol, records):
    cases = scenarios(protocol)
    expected = {s['case_id']: s for s in cases}
    observed = {}
    for record in records:
        scene = record.get('scene', {})
        key = scene.get('case_id')
        if key not in expected or scene != expected[key] or key in observed:
            raise ValueError('DUPLICATE_OR_UNREGISTERED_CONTACT_RECORD')
        if record.get('status') == 'FAILED':
            if not record.get('exception_type') or not record.get('exception_message'):
                raise ValueError('FAILED_SCENE_WITHOUT_DIAGNOSTIC')
        elif record.get('status') == 'PHYSICAL_TRAJECTORY_RECORDED_NOT_CALIBRATED':
            dt = scene['dt_seconds']
            steps = round(protocol['duration_seconds']/dt)
            trace = record.get('records', [])
            stride = max(1, round(1/(200*dt)))
            times = [i*dt for i in range(0, steps+1, stride)]
            if times[-1] != steps*dt:
                times.append(steps*dt)
            if (record.get('actual_steps') != steps or record.get('duration_seconds') != 8. or
                    [r.get('time_s') for r in trace] != times):
                raise ValueError('INCOMPLETE_CONTACT_TRAJECTORY')
            if any(type(record.get(name)) not in (int, float) or not math.isfinite(record[name])
                   for name in METRICS):
                raise ValueError('NONFINITE_CONTACT_METRIC')
        else:
            raise ValueError('UNRECOGNIZED_CONTACT_TERMINATION')
        observed[key] = record
    if set(observed) != set(expected):
        raise ValueError('MISSING_CONTACT_SCENARIOS')
    successful = [r for r in observed.values() if r['status'] != 'FAILED']
    failed = [dict(scene=r['scene'], exception_type=r['exception_type'],
                   exception_message=r['exception_message']) for r in observed.values() if r['status'] == 'FAILED']
    comparisons = []
    for left, scene in enumerate(cases):
        a = observed[scene['case_id']]
        for other in cases[left+1:]:
            same = scene['material'] == other['material'] and scene['azimuth_degrees'] == other['azimuth_degrees']
            kind = ('TIME_STEP' if same and scene['proxy'] == other['proxy'] else
                    'COLLISION_PROXY' if same and scene['dt_seconds'] == other['dt_seconds'] else None)
            if kind is None:
                continue
            b = observed[other['case_id']]
            available = a['status'] != 'FAILED' and b['status'] != 'FAILED'
            comparisons.append(dict(kind=kind, case_ids=[scene['case_id'], other['case_id']],
                status='MEASURED_DIFFERENCE' if available else 'INCOMPLETE_PAIR_FAILURE_RETAINED',
                absolute_differences={m: abs(a[m]-b[m]) for m in METRICS} if available else None))
    return dict(schema='rounded-heart-contact-full-matrix-accounting-v1',
        status='CONTACT_DATA_SUMMARIZED_NO_PHYSICAL_SUCCESS_INFERENCE',
        denominator=len(cases), recorded_scene_count=len(successful), failed_scene_count=len(failed),
        failures=failed, comparisons=comparisons,
        video_case_ids=[s['case_id'] for s in cases if s['azimuth_degrees'] == protocol['video_direction_degrees']],
        physical_validation_passed=False, calibrated_real_stopping_time_prediction=False,
        numeric_summary={m: dict(minimum=min(r[m] for r in successful),
                                 maximum=max(r[m] for r in successful)) if successful else None for m in METRICS})
