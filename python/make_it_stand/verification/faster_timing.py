"""Exclusive solver stage accounting; never infer an unmeasured paper time."""
from collections import Counter


ALGORITHM_STAGES = frozenset(('source', 'boolean', 'integration', 'derivative',
                              'intersection_motion', 'prepare_arithmetic_boundary'))
VALIDATION_STAGES = frozenset(('map_same_call_material_validation',
    'map_variable_mesh_snapshot', 'map_independent_variable_integration',
    'map_target_model_frame_identity'))
RECORDING_STAGES = frozenset(('evaluation_checkpoint', 'newton_event_record',
                              'material_attachment_serialization'))


def summarize_solver_timing(events):
    """Count exclusive durations inside the first native solve, including failures.

    Historical uncovered work stays unclassified. An interrupted prefix has no
    total solve time. These measurements do not establish convergence or an
    exactly matched reproduction of the authors' timer boundaries.
    """
    stack = []
    exclusive = Counter()
    counts = Counter()
    active = False
    completed = False
    wall = None
    for row in events:
        kind, name = row.get('kind'), row.get('stage')
        if kind == 'STAGE_STARTED':
            if not active:
                if name != 'native_newton_solver':
                    continue
                active = True
            stack.append(dict(name=name, children=0))
        elif active and kind in ('STAGE_RETURNED', 'STAGE_FAILED'):
            if not stack or stack[-1]['name'] != name:
                raise ValueError('STAGE_NESTING_MISMATCH')
            frame = stack.pop()
            duration = row.get('elapsed_ns')
            if type(duration) is not int or duration < frame['children']:
                raise ValueError('INVALID_STAGE_DURATION')
            exclusive[name] += duration-frame['children']
            counts[name] += 1
            if stack:
                stack[-1]['children'] += duration
            else:
                wall = duration
                completed = kind == 'STAGE_RETURNED'
                break
    algorithm = sum(exclusive[k] for k in ALGORITHM_STAGES)
    validation = sum(exclusive[k] for k in VALIDATION_STAGES)
    recording = sum(exclusive[k] for k in RECORDING_STAGES)
    known = algorithm+validation+recording
    return dict(schema='faster-exclusive-solver-timing-v1',
        solver_completed=completed, solver_wall_ns=wall,
        measured_algorithm_ns=algorithm, validation_ns=validation,
        recording_ns=recording,
        unclassified_ns=None if wall is None else wall-known,
        completed_stage_unclassified_ns=sum(exclusive.values())-known,
        exclusive_stage_ns=dict(exclusive), stage_counts=dict(counts),
        unfinished_stages=[frame['name'] for frame in stack],
        paper_comparable_seconds=None,
        scope='First native solver invocation; preparation and terminal seal excluded',
        limitation='Unclassified costs are not assumed to be checks or recording; '
                   'partial algorithm totals are not a time to convergence')
