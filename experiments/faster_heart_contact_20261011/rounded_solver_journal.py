"""Durable recording around the existing actual CGAL map and Faster Newton.

This numerical adapter does not admit a Heart input, validate derivatives, seal
material, or authorize simulation. Those remain separate mandatory gates.
"""
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from make_it_stand.spin import (evaluate_spin_it_faster_fullrank_quadric,
                               solve_spin_it_faster_with_external_moment_evaluator)
from make_it_stand import spin


def json_value(value):
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {'nonfinite_float': repr(value)}
    return value


def encoded(value):
    return json.dumps(json_value(value), sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def newton_backend_identity(provenance):
    """Bind complete native provenance without the Newton journal delimiter.

    The original string remains in every persisted evaluation; hashing is an
    identity encoding, not a change of backend or a geometry admission.
    """
    if not isinstance(provenance, str) or not provenance:
        raise ValueError('EMPTY_NATIVE_BACKEND_PROVENANCE')
    return 'backend-sha256-' + hashlib.sha256(provenance.encode()).hexdigest()


def save(path, value):
    with path.open('x') as stream:
        stream.write(encoded(value)+'\n'); stream.flush(); os.fsync(stream.fileno())
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class DurableFullRankMap:
    def __init__(self, request, settings, directory, *, capture_material_proof=False):
        if not isinstance(capture_material_proof, bool):
            raise ValueError('capture_material_proof must be boolean')
        self.capture_material_proof = capture_material_proof
        self.returned_mappings = []
        self.request, self.settings = request, settings
        self.directory = Path(directory); self.directory.mkdir(exist_ok=False)
        native_path = Path(spin._native_module().__file__).resolve()
        native = dict(path=str(native_path), sha256=hashlib.sha256(native_path.read_bytes()).hexdigest())
        self.source_files = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (native_path, Path(spin.__file__).resolve(), Path(__file__).resolve())}
        self.identity = digest(dict(request=asdict(request), evaluator=asdict(settings),
            source_files=self.source_files, adapter_schema=('durable-same-call-material-proof-v1'
                if capture_material_proof else 'durable-actual-fullrank-v1')))
        save(self.directory/'mapping-config.json', dict(request=asdict(request), evaluator=asdict(settings),
             native_module=native, source_files=self.source_files,
             identity=self.identity, material_acceptance='NOT_CHECKED'))
        self.mapping_stream = (self.directory/'mapping-events.jsonl').open('x')
        self.newton_stream = (self.directory/'newton-events.jsonl').open('x')
        self.count = 0

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.mapping_stream.close(); self.newton_stream.close()

    def append(self, stream, row):
        stream.write(encoded(row)+'\n'); stream.flush(); os.fsync(stream.fileno())

    def event_sink(self, event):
        self.append(self.newton_stream, dict(event))

    def __call__(self, coefficients):
        self.count += 1
        q = np.asarray(coefficients, dtype=float)
        header = dict(index=self.count, coefficients=q.tolist(), mapping_identity=self.identity)
        self.append(self.mapping_stream, {**header, 'status': 'STARTED'})
        try:
            evaluated = evaluate_spin_it_faster_fullrank_quadric(self.request, q,
                evaluator_settings=self.settings, diagnostic_stage_observer=lambda stage:
                    self.append(self.mapping_stream, {**header, 'status': 'STAGE', 'stage': stage}),
                **({'capture_material_proof': True} if self.capture_material_proof else {}))
            authority = asdict(evaluated)
            save(self.directory/f'evaluation-{self.count:06d}.json',
                 {**header, 'status': 'RETURNED_ACTUAL_MAPPING', 'evaluation': authority})
            if self.capture_material_proof:
                path = self.directory/f'evaluation-{self.count:06d}.json'
                self.returned_mappings.append(dict(index=self.count, coefficients=q.tolist(),
                    geometry_hash=evaluated.geometry_hash, path=path.name,
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            # All fields come from this one concrete evaluation. Connectivity and
            # source/branch identities are retained; they are not constant labels.
            topology = digest(dict(faces=authority['material_faces'],
                source=authority['source_scientific_hash'], branch=authority.get('branch_identity')))
            response = dict(raw_moments=evaluated.total_raw_moments,
                variable_raw_moments=evaluated.variable_raw_moments,
                derivative_defined=evaluated.derivative_defined,
                derivative_diagnostic=evaluated.derivative_diagnostic,
                geometry_hash=evaluated.geometry_hash, topology_hash=topology,
                evaluation_map_identity=self.identity,
                boolean_backend_identity=newton_backend_identity(evaluated.boolean_backend_identity),
                derivative_identity=evaluated.scientific_hash)
            if evaluated.derivative_defined:
                if evaluated.coefficient_jacobian is None:
                    raise ValueError('DEFINED_DERIVATIVE_WITHOUT_JACOBIAN')
                response['moment_jacobian'] = np.asarray(evaluated.coefficient_jacobian)
            self.append(self.mapping_stream, {**header, 'status': 'RETURNED',
                'geometry_hash': evaluated.geometry_hash, 'topology_hash': topology})
            return response
        except Exception as exc:
            self.append(self.mapping_stream, {**header, 'status': 'FAILED',
                'exception_type': type(exc).__name__, 'exception_message': str(exc)})
            raise

    def terminal_mapping(self, result):
        """Select the exact recorded final state, never the last rejected trial."""
        if not self.capture_material_proof or result.get('status') != 'CONVERGED':
            raise ValueError('TERMINAL_MAPPING_REQUIRES_CONVERGED_CAPTURE_RUN')
        q = json_value(result.get('quadric_coefficients'))
        trace = result.get('trace')
        if not trace or q is None:
            raise ValueError('TERMINAL_MAPPING_FINAL_STATE_MISSING')
        geometry = trace[-1]['material_geometry_hash']
        matches = [row for row in self.returned_mappings
                   if row['coefficients'] == q and row['geometry_hash'] == geometry]
        if not matches:
            raise ValueError('TERMINAL_MAPPING_NOT_RECORDED')
        selected = matches[-1]
        path = self.directory/selected['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != selected['sha256']:
            raise ValueError('TERMINAL_MAPPING_SAVED_CONTENT_CHANGED')
        saved = json.loads(path.read_text())
        evaluation = saved['evaluation']
        proof = json.loads(evaluation['material_proof_json'])
        if proof['geometry_hash'] != geometry or saved['mapping_identity'] != self.identity:
            raise ValueError('TERMINAL_MAPPING_PROOF_IDENTITY_MISMATCH')
        return dict(schema='same-call-terminal-material-reference-v1',
            evaluation_index=selected['index'], evaluation_file=selected['path'],
            evaluation_sha256=selected['sha256'], mapping_identity=self.identity,
            geometry_hash=geometry, material_acceptance='NOT_CHECKED',
            derivative_acceptance='NOT_CHECKED')

    def finish(self, result):
        if any(hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected
               for path, expected in self.source_files.items()):
            raise ValueError('STALE_RESULT_LOADED_MAPPING_SOURCE_CHANGED')
        rows = [json.loads(x) for x in (self.directory/'newton-events.jsonl').read_text().splitlines()]
        if rows != json_value(result['events']):
            raise ValueError('DURABLE_NEWTON_EVENT_PREFIX_MISMATCH')
        if self.capture_material_proof and result['status'] == 'CONVERGED':
            save(self.directory/'terminal-material-reference.json', self.terminal_mapping(result))
        save(self.directory/'numerical-result.json', dict(result=result,
            mapping_evaluations=self.count, mapping_identity=self.identity,
            initialization='REQUEST_SEED_ONLY_NO_EXPLICIT_QUADRICS',
            material_acceptance='NOT_CHECKED', derivative_acceptance='NOT_CHECKED'))


def run_durable_fullrank(request, evaluator_settings, directory, *, capture_material_proof=False):
    """Use the existing seeded Newton without answer-dependent initial quadrics."""
    with DurableFullRankMap(request, evaluator_settings, directory,
            capture_material_proof=capture_material_proof) as mapping:
        try:
            result = solve_spin_it_faster_with_external_moment_evaluator(request,
                mapping, event_sink=mapping.event_sink)
            mapping.finish(result)
            return result
        except Exception as exc:
            save(mapping.directory/'numerical-failure.json', dict(exception_type=type(exc).__name__,
                exception_message=str(exc), mapping_evaluations=mapping.count,
                status='INCOMPLETE_NUMERICAL_RUN_NOT_INFEASIBILITY_CERTIFICATE'))
            raise
