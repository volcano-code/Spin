"""Deterministic exclusive bundles for the strict hollow-only artifact path.

All bundle files are flat internally generated names. JSON is strict, arrays
are NPY 1.0, little endian, C-order, pickle-disabled. No overwrite or rollback.
Hashes are local integrity observations, not authenticated signatures.
"""
from __future__ import annotations

import io
import os
from pathlib import Path
import re
import stat

import numpy as np

from .spin_artifact import SealedSpinArtifact
from .spin_common import (SCOPE, byte_digest, canonical, digest, implementation_observation,
                          keys, parse, require)
from .spin_problem import MeshSnapshot, SpinArtifactProblem
from .run_evidence import _path

SCHEMA = 'sealed-spin-bundle-v1'


def _npy(array):
    stream = io.BytesIO()
    np.lib.format.write_array(stream, array, version=(1, 0), allow_pickle=False)
    return stream.getvalue()


def _name(name):
    require(isinstance(name, str) and (name in ('manifest.json', 'problem.json') or
            re.fullmatch(r'mesh_[0-9]+_(vertices|faces)\.npy', name) is not None),
            'unsupported payload name or traversal')
    return name


def _record(name, data):
    return {'name': name, 'size': len(data), 'sha256': byte_digest(data)}


class _Writer:
    def __init__(self, fd):
        self.fd = fd
        self.ordinal = 0

    def write(self, name, data):
        fd = os.open(_name(name), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                     0o600, dir_fd=self.fd)
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        return _record(name, data)

    def mesh(self, mesh):
        ordinal = self.ordinal
        self.ordinal += 1
        if mesh.explicit_empty:
            return {'explicit_empty': True}
        out = {'explicit_empty': False}
        for field, array in zip(('vertices', 'faces'), mesh.arrays()):
            record = self.write(f'mesh_{ordinal:04d}_{field}.npy', _npy(array))
            record.update({'dtype': array.dtype.str, 'shape': list(array.shape),
                           'format': 'NPY_1_0_C_ORDER'})
            out[field] = record
        return out


class _Reader:
    def __init__(self, fd):
        self.fd = fd
        self.names = set()

    def read(self, name):
        _name(name)
        require(name not in self.names, 'duplicate payload reference')
        self.names.add(name)
        fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=self.fd)
        with os.fdopen(fd, 'rb') as stream:
            before = os.fstat(stream.fileno())
            require(stat.S_ISREG(before.st_mode), 'payload must be a regular file')
            data = stream.read()
            after = os.fstat(stream.fileno())
        now = os.stat(name, dir_fd=self.fd, follow_symlinks=False)
        identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        require(identity(before) == identity(after) == identity(now), 'payload changed while reading')
        require(len(data) == before.st_size, 'partial payload read')
        return data

    def checked(self, record, *, array=False):
        fields = ('name', 'size', 'sha256', 'dtype', 'shape', 'format') if array else ('name', 'size', 'sha256')
        keys(record, fields, 'payload record')
        require(type(record['size']) is int and record['size'] >= 0, 'invalid payload size')
        require(isinstance(record['sha256'], str) and re.fullmatch('[0-9a-f]{64}', record['sha256']),
                'invalid payload hash')
        data = self.read(record['name'])
        require(_record(record['name'], data) == {k: record[k] for k in ('name', 'size', 'sha256')},
                'payload hash or size mismatch')
        return data

    def mesh(self, value):
        require(isinstance(value, dict) and type(value.get('explicit_empty')) is bool,
                'mesh needs explicit empty flag')
        if value['explicit_empty']:
            keys(value, ('explicit_empty',), 'empty mesh')
            return MeshSnapshot(None)
        keys(value, ('explicit_empty', 'vertices', 'faces'), 'mesh')
        arrays = []
        for field, dtype in (('vertices', '<f8'), ('faces', '<i4')):
            record = value[field]
            data = self.checked(record, array=True)
            require(record['dtype'] == dtype and record['format'] == 'NPY_1_0_C_ORDER', 'unsupported NPY format')
            shape = record['shape']
            require(isinstance(shape, list) and len(shape) == 2 and
                    all(type(x) is int for x in shape) and shape[0] > 0 and shape[1] == 3,
                    'invalid array shape')
            stream = io.BytesIO(data)
            require(np.lib.format.read_magic(stream) == (1, 0), 'unsupported NPY version')
            actual_shape, fortran, actual_dtype = np.lib.format.read_array_header_1_0(stream)
            require(actual_shape == tuple(shape) and not fortran and actual_dtype.str == dtype,
                    'NPY header differs from canonical array descriptor')
            require(len(data) - stream.tell() == shape[0] * 3 * np.dtype(dtype).itemsize,
                    'NPY length does not match shape')
            array = np.load(io.BytesIO(data), allow_pickle=False)
            require(_npy(array) == data, 'noncanonical or trailing NPY data')
            arrays.append(array)
        return MeshSnapshot(tuple(arrays))


def save_spin_artifact(artifact, output_directory):
    """Write a NEW directory exclusively; failures preserve every partial file.

    The manifest is written last, and no valid artifact is returned if any
    write/reload fails. Existing directories, including partial ones, are errors.
    """
    require(isinstance(artifact, SealedSpinArtifact), 'expected SealedSpinArtifact')
    root = Path(output_directory).absolute()
    _path(root.parent, directory=True)
    # Check existing directory before expensive revalidation; mkdir still provides
    # exclusive creation under races. Nothing below is ever deleted or replaced.
    if root.exists() or root.is_symlink():
        raise FileExistsError(str(root))
    artifact = artifact.revalidate()
    root.mkdir()
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        writer = _Writer(fd)
        p = artifact.problem
        def region(r):
            return {'id': r.region_id, 'density': r.density, 'mesh': writer.mesh(r.mesh)}
        problem = {'schema': 'spin-artifact-problem-payload-v1', 'identity': p.identity,
                   'original': writer.mesh(p.original),
                   'fixed_regions': [region(r) for r in p.fixed_regions],
                   'design_regions': [region(r) for r in p.design_regions],
                   'length_unit': p.length_unit, 'mass_unit': p.mass_unit,
                   'minimum_outer_wall': p.minimum_outer_wall, 'target': p.target,
                   'validation_scope': p.validation_scope, 'input_digest': p.input_digest}
        problem_record = writer.write('problem.json', canonical(problem).encode('ascii'))
        final = {'outer': writer.mesh(artifact.outer), 'material': writer.mesh(artifact.final_material),
                 'newly_carved': writer.mesh(artifact.newly_carved),
                 'design_results': [{'id': i, 'retained': writer.mesh(r), 'carved': writer.mesh(c)}
                                    for i,r,c in artifact.design_results]}
        manifest = {'schema': SCHEMA, 'validation_scope': SCOPE, 'problem': problem_record,
                    'final': final, 'solver_claim': artifact.solver_claim, 'facts': artifact.facts,
                    'scientific_digest': artifact.scientific_digest, 'evidence': artifact.evidence,
                    'evidence_digest': artifact.evidence_digest}
        require(implementation_observation() == artifact.evidence, 'stale implementation observation')
        writer.write('manifest.json', canonical(manifest).encode('ascii'))
        os.fsync(fd)
    finally:
        os.close(fd)
    # Read the actual persisted payloads, rather than blessing pre-export facts.
    return load_spin_artifact(root, expected_problem_digest=artifact.problem.input_digest)


def load_spin_artifact(bundle_directory, *, expected_problem_digest, validation_scope=SCOPE):
    """Reconstruct and revalidate, bound to the caller's original problem digest."""
    require(validation_scope == SCOPE, 'unsupported validation profile')
    require(isinstance(expected_problem_digest, str) and re.fullmatch('[0-9a-f]{64}', expected_problem_digest),
            'expected locked problem digest is required')
    root = _path(bundle_directory, directory=True)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        reader = _Reader(fd)
        manifest = parse(reader.read('manifest.json'))
        keys(manifest, ('schema', 'validation_scope', 'problem', 'final', 'solver_claim', 'facts',
                        'scientific_digest', 'evidence', 'evidence_digest'), 'manifest')
        require(manifest['schema'] == SCHEMA and manifest['validation_scope'] == SCOPE, 'unsupported bundle schema/profile')
        require(manifest['evidence'] == implementation_observation(), 'stale implementation observation')
        require(isinstance(manifest['problem'], dict) and manifest['problem'].get('name') == 'problem.json',
                'missing canonical problem payload')
        p = parse(reader.checked(manifest['problem']))
        keys(p, ('schema', 'identity', 'original', 'fixed_regions', 'design_regions', 'length_unit',
                 'mass_unit', 'minimum_outer_wall', 'target', 'validation_scope', 'input_digest'), 'problem')
        require(p['schema'] == 'spin-artifact-problem-payload-v1', 'unsupported problem schema')
        require(p['input_digest'] == expected_problem_digest, 'wrong expected problem binding')
        def regions(values):
            require(isinstance(values, list), 'material regions must be a list')
            result = []
            for r in values:
                keys(r, ('id', 'density', 'mesh'), 'material region')
                result.append((r['id'], reader.mesh(r['mesh']), r['density']))
            return result
        original = reader.mesh(p['original'])
        fixed, design = regions(p['fixed_regions']), regions(p['design_regions'])
        problem = SpinArtifactProblem(original=original, fixed_regions=fixed, design_regions=design,
            length_unit=p['length_unit'], mass_unit=p['mass_unit'], target=p['target'],
            minimum_outer_wall=p['minimum_outer_wall'], validation_scope=p['validation_scope'])
        require(problem.input_digest == expected_problem_digest and canonical(problem.identity) == canonical(p['identity']),
                'recomputed problem identity mismatch')
        final = manifest['final']
        keys(final, ('outer', 'material', 'newly_carved', 'design_results'), 'final')
        outer, material, carved = (reader.mesh(final[k]) for k in ('outer', 'material', 'newly_carved'))
        require(isinstance(final['design_results'], list), 'design results must be a list')
        results = []
        for row in final['design_results']:
            keys(row, ('id', 'retained', 'carved'), 'design result')
            results.append((row['id'], reader.mesh(row['retained']), reader.mesh(row['carved'])))
        require(set(os.listdir(fd)) == reader.names, 'unexpected or missing bundle files')
        artifact = SealedSpinArtifact.create(problem, outer=outer, design_results=results,
            final_material=material, newly_carved=carved, solver_claim=manifest['solver_claim'])
        require(canonical(artifact.facts) == canonical(manifest['facts']) and
                artifact.scientific_digest == manifest['scientific_digest'], 'recomputed scientific facts mismatch')
        require(artifact.evidence == manifest['evidence'] and artifact.evidence_digest == manifest['evidence_digest'],
                'recomputed evidence mismatch')
        return artifact
    finally:
        os.close(fd)
