"""Deterministic exclusive bundles for the strict hollow-only artifact path.

All bundle files are flat internally generated names. JSON is strict, arrays
are NPY 1.0, little endian, C-order, pickle-disabled. No overwrite or rollback.
Hashes are local integrity observations, not authenticated signatures.
"""
from __future__ import annotations

import io
from dataclasses import asdict
import os
from pathlib import Path
import re
import stat

import numpy as np

from .spin_artifact import SealedSpinArtifact
from .spin_common import (SCOPE, byte_digest, canonical, digest, implementation_observation,
                          keys, parse, require, native)
from .spin_problem import MeshSnapshot, SpinArtifactProblem
from .run_evidence import _path

SCHEMA = 'sealed-spin-bundle-v2'
LEGACY_SCHEMA = 'sealed-spin-bundle-v1'
OWNED_SCHEMA = 'sealed-spin-bundle-v3-owned-operations'


def _owned_mesh_operands(mesh):
    if mesh.explicit_empty:
        return ()
    region = mesh.exact_region()
    operands = getattr(region, 'construction_operands', ())
    if region.facts['algorithm'] == 'strict-interior-nef-difference-shared-boundary-v1':
        require(len(operands) == 2, 'owned operation replay metadata unavailable in this native build')
        return tuple(operands)
    require(not operands, 'unsupported owned region operation')
    return ()


def _npy(array):
    stream = io.BytesIO()
    np.lib.format.write_array(stream, array, version=(1, 0), allow_pickle=False)
    return stream.getvalue()


def _name(name):
    require(isinstance(name, str) and (name in ('manifest.json', 'problem.json') or
            re.fullmatch(r'mesh_[0-9]+_(vertices|faces)\.npy', name) is not None or
            re.fullmatch(r'mesh_[0-9]+_rational_vertices\.json', name) is not None),
            'unsupported payload name or traversal')
    return name


def _record(name, data):
    return {'name': name, 'size': len(data), 'sha256': byte_digest(data)}


class _Writer:
    def __init__(self, fd, *, schema=SCHEMA):
        self.fd = fd
        self.ordinal = 0
        self.schema = schema

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
            out = {'explicit_empty': True, 'rational_vertices': None}
            if self.schema == OWNED_SCHEMA:
                out['construction'] = None
            return out
        out = {'explicit_empty': False, 'rational_vertices': None}
        for field, array in zip(('vertices', 'faces'), mesh.arrays()):
            record = self.write(f'mesh_{ordinal:04d}_{field}.npy', _npy(array))
            record.update({'dtype': array.dtype.str, 'shape': list(array.shape),
                           'format': 'NPY_1_0_C_ORDER'})
            out[field] = record
        exact = mesh.rational_vertices
        if exact is not None:
            out['rational_vertices'] = self.write(
                f'mesh_{ordinal:04d}_rational_vertices.json',
                canonical([list(row) for row in exact]).encode('ascii'))
        operands = _owned_mesh_operands(mesh)
        require(not operands or self.schema == OWNED_SCHEMA, 'owned operations require the v3 bundle schema')
        if self.schema == OWNED_SCHEMA:
            out['construction'] = (dict(operation='STRICT_INTERIOR_DIFFERENCE_V1',
                operands=[self.mesh(MeshSnapshot.from_exact_region(region)) for region in operands])
                if operands else None)
        return out


class _Reader:
    def __init__(self, fd):
        self.fd = fd
        self.names = set()
        # Identical exact meshes can be referenced in multiple bundle roles.
        # Cache immutable snapshots by the full canonical payload identity;
        # every reference is still independently read and integrity-checked.
        self._exact_mesh_cache = {}

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

    def mesh(self, value, *, schema):
        require(isinstance(value, dict) and type(value.get('explicit_empty')) is bool,
                'mesh needs explicit empty flag')
        v2 = schema in (SCHEMA, OWNED_SCHEMA)
        owned = schema == OWNED_SCHEMA
        expected = ('explicit_empty', 'rational_vertices') if v2 else ('explicit_empty',)
        if owned:
            expected += ('construction',)
        if value['explicit_empty']:
            keys(value, expected, 'empty mesh')
            require(not v2 or value['rational_vertices'] is None,
                    'empty mesh cannot carry rational vertices')
            require(not owned or value['construction'] is None, 'empty mesh cannot carry an operation')
            return MeshSnapshot(None)
        fields = ('explicit_empty', 'vertices', 'faces', 'rational_vertices') if v2 else (
            'explicit_empty', 'vertices', 'faces')
        if owned:
            fields += ('construction',)
        keys(value, fields, 'mesh')
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
        exact_data, exact_vertices = None, None
        if v2 and value['rational_vertices'] is not None:
            exact_data = self.checked(value['rational_vertices'])
            exact_vertices = parse(exact_data.decode('ascii'))
            require(canonical(exact_vertices).encode('ascii') == exact_data,
                    'noncanonical rational vertex payload')
        recipe = value.get('construction') if owned else None
        operands, operation_identity = None, None
        if recipe is not None:
            keys(recipe, ('operation', 'operands'), 'owned mesh construction')
            require(recipe['operation'] == 'STRICT_INTERIOR_DIFFERENCE_V1', 'unsupported owned mesh operation')
            require(type(recipe['operands']) is list and len(recipe['operands']) == 2,
                    'owned difference requires exactly two actual operands')
            operands = [self.mesh(item, schema=schema) for item in recipe['operands']]
            operation_identity = (recipe['operation'], tuple(item.exact_region().facts['region_digest'] for item in operands))
        cache_key = (arrays[0].tobytes(order='C'), arrays[1].tobytes(order='C'), exact_data, operation_identity)
        cached = self._exact_mesh_cache.get(cache_key)
        if cached is not None:
            return cached
        if operands is not None:
            region = native().make_exact_interior_difference_region(*(item.exact_region() for item in operands))
            exact = MeshSnapshot.from_exact_region(region)
        elif exact_vertices is None:
            exact = MeshSnapshot(tuple(arrays))
        else:
            exact = MeshSnapshot.from_rationals(exact_vertices, arrays[1])
        require(all(a.tobytes(order='C') == b.tobytes(order='C') for a,b in zip(exact.arrays(), arrays)),
                'rational payload or owned operation does not match derived binary64 vertices/faces')
        actual_exact = (None if exact.rational_vertices is None else
                        canonical([list(row) for row in exact.rational_vertices]).encode('ascii'))
        require(actual_exact == exact_data, 'owned operation does not match the actual rational payload')
        self._exact_mesh_cache[cache_key] = exact
        return exact


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
    from .spin_artifact import _reuse_created_artifact_validation
    if not _reuse_created_artifact_validation(artifact):
        artifact = artifact.revalidate()
    root.mkdir()
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        p = artifact.problem
        snapshots = [p.original, artifact.outer, artifact.final_material, artifact.newly_carved]
        snapshots += [r.mesh for r in p.fixed_regions+p.design_regions]
        snapshots += [mesh for _, retained, carved in artifact.design_results for mesh in (retained, carved)]
        schema = OWNED_SCHEMA if any(_owned_mesh_operands(mesh) for mesh in snapshots) else SCHEMA
        writer = _Writer(fd, schema=schema)
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
        manifest = {'schema': schema, 'validation_scope': SCOPE, 'problem': problem_record,
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
    artifact, _ = _read_and_validate_spin_artifact(bundle_directory,
        expected_problem_digest=expected_problem_digest, validation_scope=validation_scope)
    return artifact


STATIC_OUTPUT_SCHEMA = 'spin2014-static-material-attachment-bundle-v1'
ROTATIONAL_OUTPUT_SCHEMA = 'spin2014-rotational-material-attachment-bundle-v1'


def _save_spin2014_authority_point_attachment(point, settings, output_directory, *,
                                             expected_request_identity, expected_context_identity,
                                             arithmetic_backend='FRACTION_V1'):
    """Persist an owned authority point using existing bundle machinery.

    This internal attachment accepts diagnostic points and NEVER grants a final
    optimization/manufacturing seal. Final run admission remains a separate gate.
    Prepare once here; then read files and independently prepare again below.
    """
    from ._spin2014_authority import prepare_authority_material_boundaries, authority_arithmetic_identity
    arithmetic=authority_arithmetic_identity(arithmetic_backend)
    root = _rotational_output_destination(output_directory)
    # Canonical snapshot prevents caller mutation between preparation and writing.
    accelerated=arithmetic_backend!='FRACTION_V1'
    source_value={'schema':'spin2014-authority-point-source-v2' if accelerated else 'spin2014-authority-point-source-v1',
                  'point':point,'settings':settings}
    if accelerated:
        source_value['arithmetic_identity']=arithmetic
    source = parse(canonical(source_value))
    regions,audit = prepare_authority_material_boundaries(source['point'],source['settings'],
        expected_request_identity=expected_request_identity,
        expected_context_identity=expected_context_identity,arithmetic_backend=arithmetic_backend)
    root.mkdir()
    fd = os.open(root,os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        writer = _Writer(fd)
        record = writer.write('problem.json',canonical(source).encode('ascii'))
        meshes = {name:writer.mesh(mesh) for name,mesh in regions.items()}
        writer.write('manifest.json',canonical({
            'schema':'spin2014-authority-point-attachment-v2' if accelerated else 'spin2014-authority-point-attachment-v1','source':record,
            'meshes':meshes,'audit':audit}).encode('ascii'))
        os.fsync(fd)
    finally:
        os.close(fd)
    return _load_spin2014_authority_point_attachment(root,
        expected_point_digest=audit['source_point_digest'],
        expected_request_identity=expected_request_identity,
        expected_context_identity=expected_context_identity)


def _load_spin2014_authority_point_attachment(bundle_directory, *, expected_point_digest,
                                             expected_request_identity, expected_context_identity):
    """Independent file-backed reconstruction, not a cached preparation result."""
    from ._spin2014_authority import prepare_authority_material_boundaries, authority_arithmetic_identity
    require(isinstance(expected_point_digest,str) and
            re.fullmatch('[0-9a-f]{64}',expected_point_digest),'expected authority point digest required')
    root = _path(bundle_directory,directory=True)
    fd = os.open(root,os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        reader = _Reader(fd)
        manifest = parse(reader.read('manifest.json'))
        keys(manifest,('schema','source','meshes','audit'),'authority point manifest')
        require(manifest['schema'] in ('spin2014-authority-point-attachment-v1','spin2014-authority-point-attachment-v2'),
                'unsupported authority point attachment')
        require(manifest['source'].get('name')=='problem.json','authority source payload missing')
        source = parse(reader.checked(manifest['source']))
        accelerated=manifest['schema']=='spin2014-authority-point-attachment-v2'
        keys(source,('schema','point','settings','arithmetic_identity') if accelerated else ('schema','point','settings'),'authority point source')
        arithmetic_backend='GMPY2_2_2_1_EXACT_V1' if accelerated else 'FRACTION_V1'
        if accelerated:
            require(canonical(source['arithmetic_identity'])==canonical(authority_arithmetic_identity(arithmetic_backend)),
                    'authority attachment arithmetic identity mismatch')
        require(source['schema']==('spin2014-authority-point-source-v2' if accelerated else 'spin2014-authority-point-source-v1') and
                digest(source['point'])==expected_point_digest,'authority source binding mismatch')
        names = ('outer','fixed','design','retained','cavity','material')
        keys(manifest['meshes'],names,'authority point meshes')
        meshes = {name:reader.mesh(manifest['meshes'][name],schema=SCHEMA) for name in names}
        require(set(os.listdir(fd))==reader.names,'unexpected or missing authority attachment files')
        regions,audit = prepare_authority_material_boundaries(source['point'],source['settings'],
            expected_request_identity=expected_request_identity,
            expected_context_identity=expected_context_identity,arithmetic_backend=arithmetic_backend)
        for name in names:
            require(canonical(meshes[name].identity)==canonical(regions[name].identity),
                    'saved authority boundary differs from reconstructed source: '+name)
        require(canonical(audit)==canonical(manifest['audit']),
                'saved authority audit differs from independent reconstruction')
        return regions,audit
    finally:
        os.close(fd)


def _rotational_source_run(value):
    """Restore the complete saved run, including failed trials, without solving."""
    from ..spin import (Spin2014QuadricMaterialRun, Spin2014QuadricObjectiveConfig,
                        SpinItFasterRequest, SpinSolverConfig, _fullrank_solve_result_from_native)
    from .faster_material import _json_value
    keys(value, ('objective_profile', 'objective_identity', 'objective_settings',
                 'measured_request', 'fixed_material_density', 'solve', 'material_export_complete',
                 'material_export_diagnostic', 'cavity_vertices', 'cavity_faces',
                 'cavity_raw_moments'), 'rotational source run')

    def restore(item):
        if isinstance(item, dict):
            if '$nonfinite_float' in item:
                keys(item, ('$nonfinite_float',), 'nonfinite diagnostic')
                require(item['$nonfinite_float'] in ('nan', '+inf', '-inf'),
                        'invalid nonfinite diagnostic encoding')
                return float({'nan':'nan', '+inf':'inf', '-inf':'-inf'}[item['$nonfinite_float']])
            return {key:restore(val) for key,val in item.items()}
        if isinstance(item, list):
            return [restore(val) for val in item]
        return item

    request = dict(value['measured_request'])
    request['settings'] = SpinSolverConfig(**request['settings'])
    out = Spin2014QuadricMaterialRun(
        objective_profile=value['objective_profile'], objective_identity=value['objective_identity'],
        objective_settings=Spin2014QuadricObjectiveConfig(**value['objective_settings']),
        measured_request=SpinItFasterRequest(**request),
        fixed_material_density=value['fixed_material_density'],
        solve=_fullrank_solve_result_from_native(restore(value['solve'])),
        material_export_complete=value['material_export_complete'],
        material_export_diagnostic=value['material_export_diagnostic'],
        cavity_vertices=value['cavity_vertices'], cavity_faces=value['cavity_faces'],
        cavity_raw_moments=value['cavity_raw_moments'])
    require(canonical(_json_value(asdict(out))) == canonical(value),
            'rotational source parser changed or omitted saved evidence')
    return out


def _rotational_output_destination(output_directory):
    root = Path(output_directory).absolute()
    _path(root.parent, directory=True)
    if root.exists() or root.is_symlink():
        raise FileExistsError(str(root))
    return root


def prepare_and_save_spin2014_quadric_material_output(source_run, output_directory, *, evaluator_settings):
    """Prepare once, persist that material, then independently reload/audit.

    Unlike the legacy prepared-output entry, the caller supplies no material
    object or approval flag. The saved source and all six boundaries are
    independently reconstructed by the unchanged post-save loader.
    """
    from .faster_material import prepare_spin2014_quadric_material_output
    root = _rotational_output_destination(output_directory)
    current = prepare_spin2014_quadric_material_output(source_run, evaluator_settings=evaluator_settings)
    return _write_prepared_spin2014_quadric_material_output(
        current, source_run, root, evaluator_settings=evaluator_settings)


def save_spin2014_quadric_material_output(output, source_run, output_directory, *, evaluator_settings):
    """Exclusively persist the actual rotational material, then reload/audit it.

    Uses the existing exact mesh bundle I/O, not the manufacturing seal schema.
    The complete numerical source and explicit evaluator are saved, including
    failures and timings. No optimization, repair or backend fallback occurs.
    """
    from .faster_material import (Spin2014QuadricMaterialOutput, _json_value,
                                 prepare_spin2014_quadric_material_output)
    require(isinstance(output, Spin2014QuadricMaterialOutput), 'expected rotational material output')
    root = _rotational_output_destination(output_directory)
    source = _json_value(asdict(source_run))
    require(digest(source) == output.identity['source_run_sha256'], 'rotational source binding mismatch')
    current = prepare_spin2014_quadric_material_output(source_run, evaluator_settings=evaluator_settings)
    require(canonical(current.identity) == canonical(output.identity) and
            current.scientific_digest == output.scientific_digest,
            'rotational output identity or audit differs from its actual source')
    return _write_prepared_spin2014_quadric_material_output(
        current, source_run, root, evaluator_settings=evaluator_settings)


def _write_prepared_spin2014_quadric_material_output(current, source_run, root, *, evaluator_settings):
    """Internal writer; public entries have just prepared this exact source."""
    from .faster_material import _json_value
    source = _json_value(asdict(source_run))
    require(digest(source) == current.identity['source_run_sha256'], 'rotational source binding mismatch')
    root.mkdir()
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        writer = _Writer(fd)
        source_record = writer.write('problem.json', canonical({
            'schema':'spin2014-rotational-material-source-v1', 'run':source,
            'evaluator_settings':_json_value(evaluator_settings.identity_dict())}).encode('ascii'))
        meshes = {name:writer.mesh(getattr(current,name)) for name in
                  ('outer', 'fixed', 'design', 'retained', 'cavity', 'material')}
        writer.write('manifest.json', canonical({
            'schema':ROTATIONAL_OUTPUT_SCHEMA, 'source':source_record, 'meshes':meshes,
            'identity':current.identity, 'scientific_digest':current.scientific_digest}).encode('ascii'))
        os.fsync(fd)
    finally:
        os.close(fd)
    return load_spin2014_quadric_material_output(root,
        expected_source_run_digest=current.identity['source_run_sha256'])


def load_spin2014_quadric_material_output(bundle_directory, *, expected_source_run_digest):
    """Read/authenticate all saved meshes and independently recompute the audit.

    The saved source is reconstructed but never optimized. Exact coordinates
    remain authoritative; a failed floating derivation is not repaired or
    promoted to manufacturing. Caller-supplied audit booleans cannot pass.
    """
    from ..spin import SpinEvaluatorConfig
    from .faster_material import prepare_spin2014_quadric_material_output
    require(isinstance(expected_source_run_digest,str) and
            re.fullmatch('[0-9a-f]{64}', expected_source_run_digest), 'expected source binding is required')
    root = _path(bundle_directory,directory=True)
    fd = os.open(root,os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        reader = _Reader(fd)
        manifest = parse(reader.read('manifest.json'))
        keys(manifest, ('schema','source','meshes','identity','scientific_digest'), 'rotational material manifest')
        require(manifest['schema'] == ROTATIONAL_OUTPUT_SCHEMA, 'unsupported rotational attachment schema')
        require(manifest['source'].get('name') == 'problem.json', 'missing rotational source payload')
        source = parse(reader.checked(manifest['source']))
        keys(source, ('schema','run','evaluator_settings'), 'rotational material source')
        require(source['schema'] == 'spin2014-rotational-material-source-v1' and
                digest(source['run']) == expected_source_run_digest, 'wrong expected source run binding')
        names = ('outer','fixed','design','retained','cavity','material')
        keys(manifest['meshes'],names,'rotational material meshes')
        meshes = {name:reader.mesh(manifest['meshes'][name],schema=SCHEMA) for name in names}
        require(set(os.listdir(fd)) == reader.names, 'unexpected or missing rotational material payload files')
        current = prepare_spin2014_quadric_material_output(_rotational_source_run(source['run']),
            evaluator_settings=SpinEvaluatorConfig(**source['evaluator_settings']))
        for name,mesh in meshes.items():
            if canonical(mesh.identity) != canonical(getattr(current,name).identity):
                def geometry(snapshot):
                    vertices,faces=snapshot.arrays()
                    return {'vertices':vertices.tolist(),'faces':faces.tolist(),
                            'rational_vertices':snapshot.rational_vertices}
                raise ValueError('saved rotational material differs from actual source: '+name+
                    '; boundary witness '+canonical({'saved':geometry(mesh),
                        'reconstructed':geometry(getattr(current,name))}))
        require(canonical(current.identity) == canonical(manifest['identity']) and
                current.scientific_digest == manifest['scientific_digest'],
                'recomputed rotational material identity or audit mismatch')
        return current
    finally:
        os.close(fd)


def _static_source_run(value):
    """Restore a saved attachment only; never initialize or call a solver."""
    from ..spin import (Spin2014StaticQuadricMaterialRun, Spin2014QuadricObjectiveConfig,
                        SpinItFasterRequest, SpinSolverConfig)
    keys(value, ('objective_profile', 'objective_identity', 'objective_settings',
                 'measured_request', 'fixed_material_density', 'target_length',
                 'initialization', 'initial_coefficients', 'candidate', 'squared_height',
                 'material_verified'), 'static source run')
    require(value['material_verified'] is False, 'saved numerical run cannot claim material approval')
    request = dict(value['measured_request'])
    request['settings'] = SpinSolverConfig(**request['settings'])
    return Spin2014StaticQuadricMaterialRun(
        objective_profile=value['objective_profile'], objective_identity=value['objective_identity'],
        objective_settings=Spin2014QuadricObjectiveConfig(**value['objective_settings']),
        measured_request=SpinItFasterRequest(**request),
        fixed_material_density=value['fixed_material_density'], target_length=value['target_length'],
        initialization=value['initialization'], initial_coefficients=tuple(value['initial_coefficients']),
        candidate=value['candidate'], squared_height=value['squared_height'])


def save_spin2014_static_material_output(output, source_run, output_directory):
    """Persist a versioned geometry attachment using the existing bundle I/O.

    This is deliberately not the sealed bundle schema: wall, support, fresh
    solver evidence and independent repeat gates remain separate requirements.
    Saved files are reloaded and audited before returning any output object.
    """
    from .faster_material import (Spin2014StaticMaterialOutput, _json_value,
                                 prepare_spin2014_static_material_output)
    require(isinstance(output, Spin2014StaticMaterialOutput), 'expected a static material output')
    root = Path(output_directory).absolute()
    _path(root.parent, directory=True)
    if root.exists() or root.is_symlink():
        raise FileExistsError(str(root))
    source = _json_value(asdict(source_run))
    require(digest(source) == output.identity['source_run_sha256'], 'static output source run binding mismatch')
    current = prepare_spin2014_static_material_output(source_run)
    require(current.scientific_digest == output.scientific_digest,
            'static output identity does not match its reaudited source run')
    root.mkdir()
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        writer = _Writer(fd)
        source_record = writer.write('problem.json', canonical({
            'schema': 'spin2014-static-material-source-v1', 'run': source}).encode('ascii'))
        meshes = {name: writer.mesh(getattr(current, name)) for name in
                  ('outer', 'fixed', 'design', 'retained', 'cavity', 'material')}
        writer.write('manifest.json', canonical({
            'schema': STATIC_OUTPUT_SCHEMA, 'source': source_record, 'meshes': meshes,
            'identity': current.identity, 'scientific_digest': current.scientific_digest}).encode('ascii'))
        os.fsync(fd)
    finally:
        os.close(fd)
    return load_spin2014_static_material_output(root,
        expected_source_run_digest=output.identity['source_run_sha256'])


def load_spin2014_static_material_output(bundle_directory, *, expected_source_run_digest):
    """Recheck actual saved geometry, exact partitions and independent moments.

    Source geometry is regenerated for binding, not optimized. Every persisted
    mesh is independently read/authenticated and must equal that actual source
    construction. A caller-provided pass flag or a modified audit cannot grant
    material, support, manufacturing or functional completion.
    """
    from .faster_material import prepare_spin2014_static_material_output
    require(isinstance(expected_source_run_digest, str) and
            re.fullmatch('[0-9a-f]{64}', expected_source_run_digest), 'expected source binding is required')
    root = _path(bundle_directory, directory=True)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        reader = _Reader(fd)
        manifest = parse(reader.read('manifest.json'))
        keys(manifest, ('schema', 'source', 'meshes', 'identity', 'scientific_digest'), 'static material manifest')
        require(manifest['schema'] == STATIC_OUTPUT_SCHEMA, 'unsupported static material attachment schema')
        require(manifest['source'].get('name') == 'problem.json', 'missing static source payload')
        source = parse(reader.checked(manifest['source']))
        keys(source, ('schema', 'run'), 'static material source')
        require(source['schema'] == 'spin2014-static-material-source-v1' and
                digest(source['run']) == expected_source_run_digest, 'wrong expected source run binding')
        names = ('outer', 'fixed', 'design', 'retained', 'cavity', 'material')
        keys(manifest['meshes'], names, 'static material meshes')
        meshes = {name: reader.mesh(manifest['meshes'][name], schema=SCHEMA) for name in names}
        require(set(os.listdir(fd)) == reader.names, 'unexpected or missing static material payload files')
        current = prepare_spin2014_static_material_output(_static_source_run(source['run']))
        for name, mesh in meshes.items():
            if canonical(mesh.identity) != canonical(getattr(current, name).identity):
                # Retain actual reloaded/reconstructed boundaries in the
                # diagnostic. Do not replace either or accept mere equality
                # of moments as proof that the entity is the same.
                def geometry(snapshot):
                    v, f = snapshot.arrays()
                    return {'vertices': v.tolist(), 'faces': f.tolist(),
                            'rational_vertices': snapshot.rational_vertices}
                raise ValueError('saved static material mesh differs from its actual source: ' + name+
                    '; boundary witness '+canonical({'saved':geometry(mesh),
                        'reconstructed':geometry(getattr(current,name))}))
        require(canonical(current.identity) == canonical(manifest['identity']) and
                current.scientific_digest == manifest['scientific_digest'],
                'recomputed static material identity or audit mismatch')
        return current
    finally:
        os.close(fd)


def reaudit_spin_artifact(bundle_directory, *, expected_problem_digest,
                         expected_scientific_digest, validation_scope=SCOPE):
    """Recheck saved geometry using current validators, without editing history.

    Requires exactly identical recomputed scientific facts, not relaxed
    tolerances. The returned object has a NEW checker observation; the receipt
    retains the old observation and file identity. It does not refresh the
    historical solver run, attest its build, or add manufacturing approval.
    ``load_spin_artifact`` remains strict about checker identity.
    """
    require(isinstance(expected_scientific_digest, str) and
            re.fullmatch('[0-9a-f]{64}', expected_scientific_digest),
            'expected scientific digest binding is required')
    return _read_and_validate_spin_artifact(bundle_directory,
        expected_problem_digest=expected_problem_digest, validation_scope=validation_scope,
        reaudit_scientific_digest=expected_scientific_digest)


def _read_and_validate_spin_artifact(bundle_directory, *, expected_problem_digest,
                                     validation_scope, reaudit_scientific_digest=None):
    require(validation_scope == SCOPE, 'unsupported validation profile')
    require(isinstance(expected_problem_digest, str) and re.fullmatch('[0-9a-f]{64}', expected_problem_digest),
            'expected locked problem digest is required')
    root = _path(bundle_directory, directory=True)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        reader = _Reader(fd)
        manifest_bytes = reader.read('manifest.json')
        manifest = parse(manifest_bytes)
        keys(manifest, ('schema', 'validation_scope', 'problem', 'final', 'solver_claim', 'facts',
                        'scientific_digest', 'evidence', 'evidence_digest'), 'manifest')
        require(manifest['schema'] in (SCHEMA, LEGACY_SCHEMA, OWNED_SCHEMA) and manifest['validation_scope'] == SCOPE, 'unsupported bundle schema/profile')
        require(manifest['evidence_digest'] == digest({
            'scientific_digest': manifest['scientific_digest'],
            'observation': manifest['evidence']}), 'source evidence digest mismatch')
        if reaudit_scientific_digest is None:
            require(manifest['evidence'] == implementation_observation(), 'stale implementation observation')
        else:
            require(manifest['scientific_digest'] == reaudit_scientific_digest,
                    'wrong expected scientific digest binding')
        require(isinstance(manifest['problem'], dict) and manifest['problem'].get('name') == 'problem.json',
                'missing canonical problem payload')
        schema = manifest['schema']
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
                result.append((r['id'], reader.mesh(r['mesh'], schema=schema), r['density']))
            return result
        original = reader.mesh(p['original'], schema=schema)
        fixed, design = regions(p['fixed_regions']), regions(p['design_regions'])
        problem = SpinArtifactProblem(original=original, fixed_regions=fixed, design_regions=design,
            length_unit=p['length_unit'], mass_unit=p['mass_unit'], target=p['target'],
            minimum_outer_wall=p['minimum_outer_wall'], validation_scope=p['validation_scope'])
        require(problem.input_digest == expected_problem_digest and canonical(problem.identity) == canonical(p['identity']),
                'recomputed problem identity mismatch')
        final = manifest['final']
        keys(final, ('outer', 'material', 'newly_carved', 'design_results'), 'final')
        outer, material, carved = (reader.mesh(final[k], schema=schema) for k in ('outer', 'material', 'newly_carved'))
        require(isinstance(final['design_results'], list), 'design results must be a list')
        results = []
        for row in final['design_results']:
            keys(row, ('id', 'retained', 'carved'), 'design result')
            results.append((row['id'], reader.mesh(row['retained'], schema=schema), reader.mesh(row['carved'], schema=schema)))
        require(set(os.listdir(fd)) == reader.names, 'unexpected or missing bundle files')
        artifact = SealedSpinArtifact.create(problem, outer=outer, design_results=results,
            final_material=material, newly_carved=carved, solver_claim=manifest['solver_claim'],
            _reuse_problem_validation=True)
        require(canonical(artifact.facts) == canonical(manifest['facts']) and
                artifact.scientific_digest == manifest['scientific_digest'], 'recomputed scientific facts mismatch')
        if reaudit_scientific_digest is None:
            require(artifact.evidence == manifest['evidence'] and artifact.evidence_digest == manifest['evidence_digest'],
                    'recomputed evidence mismatch')
            return artifact, None
        receipt = {'schema': 'sealed-spin-material-reaudit-v1',
                   'status': 'MATERIAL_REAUDIT_IDENTICAL_FACTS',
                   'source_manifest_sha256': byte_digest(manifest_bytes),
                   'problem_digest': artifact.problem.input_digest,
                   'scientific_digest': artifact.scientific_digest,
                   'source_evidence_digest': manifest['evidence_digest'],
                   'source_observation': manifest['evidence'],
                   'current_evidence_digest': artifact.evidence_digest,
                   'current_observation': artifact.evidence,
                   'solver_freshness': 'NOT_REASSESSED',
                   'manufacturing_approval': 'NOT_GRANTED'}
        return artifact, receipt
    finally:
        os.close(fd)
