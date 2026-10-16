"""Diagnostic bridge decorator: persist literal operands before native entry.

The caller supplies a durable publisher (RunJournal.publish). This module does
not change the Boolean, classify a timeout, or provide a fallback backend.
"""
import copy
import hashlib
import json


def replay_event_projection(value):
    from make_it_stand.verification.cherchi_inprocess import _json_value
    from make_it_stand.verification.cherchi_run_evidence import scientific_projection
    return scientific_projection(json.loads(json.dumps(value, default=_json_value)))


class StopBeforeBoolean(BaseException):
    pass


class CapturingBridge:
    def __init__(self, bridge, publish, context, *, stop_at=None):
        if not callable(publish) or not callable(context):
            raise TypeError('capture needs publisher and context provider')
        if stop_at is not None and (type(stop_at) is not int or stop_at < 1):
            raise ValueError('invalid stop index')
        self.bridge, self.publish, self.context = bridge, publish, context
        self.stop_at = stop_at
        self.calls = 0

    def __getattr__(self, name):
        return getattr(self.bridge, name)

    def boolean(self, design_vertices, design_faces, source_vertices, source_faces, operation):
        self.calls += 1
        body = copy.deepcopy(dict(schema='cherchi-pre-boolean-operands-v1',
            boolean_call_index=self.calls, context=self.context(),
            design_vertices=design_vertices, design_faces=design_faces,
            source_vertices=source_vertices, source_faces=source_faces,
            operation=operation))
        encoded = json.dumps(body, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
        body['operand_packet_sha256'] = hashlib.sha256(encoded).hexdigest()
        self.publish(f'boolean-input-{self.calls:04d}.json', body)
        if self.calls == self.stop_at:
            raise StopBeforeBoolean(f'DIAGNOSTIC_STOP_BEFORE_BOOLEAN:{self.calls}')
        return self.bridge.boolean(design_vertices, design_faces, source_vertices, source_faces, operation)
