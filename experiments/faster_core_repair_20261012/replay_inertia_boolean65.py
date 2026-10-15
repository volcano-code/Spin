"""Seed replay with unchanged solver/backend, stopping before evaluation65 Boolean."""
from dataclasses import asdict
import json
from pathlib import Path
import sys
import numpy as np
import os

from make_it_stand import spin
from make_it_stand.verification import cherchi_inprocess as ci
from make_it_stand.verification.cherchi_run_journal import RunJournal
from make_it_stand.verification.cherchi_run_evidence import scientific_projection
from boolean_operand_capture import CapturingBridge, StopBeforeBoolean
from run_core_benchmark import main as benchmark


def main(output, original):
    original_run = original/'output/run'
    expected_events = [scientific_projection(x['event']) for x in
        map(json.loads, (original_run/'events.jsonl').read_text().splitlines())
        if x.get('kind') == 'NEWTON_EVENT']
    output.mkdir(parents=True, exist_ok=False)
    original_publish, original_emit = RunJournal.publish, RunJournal.emit
    original_call = ci.InProcessCherchiEvaluator.__call__
    state = dict(events=0, materials=0, evaluations=0)
    (output/'capture').mkdir()
    with RunJournal(output/'capture', json_value=ci._json_value) as journal:
        def publish(self, name, value):
            if name.startswith('evaluation-'):
                expected = json.loads((original_run/name).read_text())
                if json.loads(json.dumps(value, default=ci._json_value)) != expected:
                    raise ValueError('REPLAY_MATERIAL_PREFIX_MISMATCH:'+name)
                state['materials'] += 1
            return original_publish(self, name, value)
        def emit(self, row):
            if row.get('kind') == 'NEWTON_EVENT':
                index = state['events']
                if index >= len(expected_events) or scientific_projection(row['event']) != expected_events[index]:
                    raise ValueError('REPLAY_SCIENTIFIC_EVENT_PREFIX_MISMATCH:'+str(index))
                state['events'] += 1
            return original_emit(self, row)
        def call(self, q):
            state['evaluations'] += 1
            context = dict(evaluation_index=state['evaluations'], coefficients=list(q),
                request=asdict(self.request), backend_identity=self.backend_identity,
                arithmetic_identity=self.arithmetic_identity, context_identity=self.context_identity)
            if state['evaluations'] == 65:
                base = next(e for e in reversed(expected_events) if e['type'] == 'newton_linear_solve')
                before = np.asarray(base['quadric_coefficients'], dtype=float).copy()
                after = np.asarray(q, dtype=float).copy()
                before[7:] *= 2.; after[7:] *= 2.
                before /= np.linalg.norm(before); after /= np.linalg.norm(after)
                dot = float(before @ after)
                context['geodesic_step'] = float(np.arctan2(np.linalg.norm(after-dot*before), dot))
                context['step_evidence'] = 'measured dual-sphere arc from saved Newton base to actual trial; not native event'
                context['base_coefficients'] = base['quadric_coefficients']
            bridge = self.bridge
            original_stage = self._stage
            def stage(name, operation):
                value = original_stage(name, operation)
                if name == 'source':
                    context['source_motion'] = json.loads(json.dumps(value, default=ci._json_value))
                return value
            def persist(name, body):
                if state['evaluations'] != 65:
                    return
                if state['events'] != len(expected_events):
                    raise ValueError('REPLAY_EVENT_PREFIX_INCOMPLETE')
                expected_count = len(list(original_run.glob('evaluation-*.json')))
                if state['materials'] != expected_count:
                    raise ValueError('REPLAY_MATERIAL_PREFIX_INCOMPLETE')
                journal.publish('boolean-evaluation65.json', body)
                journal.publish('prefix-comparison.json', dict(state, passed=True))
                journal.publish('report.json', dict(status='CAPTURED_BEFORE_BOOLEAN', prefix_verified=True,
                    original_job=str(original), optimization_complete=False))
                journal.close()
                # pybind wraps even Python BaseException in std::exception;
                # Newton handles that as a rejected trial. End this isolated
                # diagnostic worker only after durable capture has closed.
                os._exit(0)
            self.bridge = CapturingBridge(bridge, persist, lambda: context)
            self._stage = stage
            try:
                return original_call(self, q)
            finally:
                self.bridge = bridge
                self._stage = original_stage
        RunJournal.publish, RunJournal.emit = publish, emit
        ci.InProcessCherchiEvaluator.__call__ = call
        try:
            benchmark(output/'replay', 'small', 1, 'NATIVE_GMP_SURFACE_V1',
                'EXACT_RATIONAL_FIXED_CLASSIFIED_V1', 'EVALUATION_CHECKPOINT_V1',
                'LAYERED_SAME_CALL_V1', 'CONSTRUCTOR_PROOF_REUSE_V1', 'OWNED_INTERIOR_DIFFERENCE_V1',
                spin.Spin2014QuadricObjectiveConfig(objective='TOP_INERTIA_ONLY'))
        except StopBeforeBoolean:
            journal.publish('report.json', dict(status='CAPTURED_BEFORE_BOOLEAN', prefix_verified=True,
                original_job=str(original), optimization_complete=False))
        finally:
            RunJournal.publish, RunJournal.emit = original_publish, original_emit
            ci.InProcessCherchiEvaluator.__call__ = original_call


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
