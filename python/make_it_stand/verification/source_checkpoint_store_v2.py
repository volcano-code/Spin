"""Opt-in v2 identity manifest for the append-only stage carrier.

Retains the SPCKPT01 framing and 1 MiB header bound, with a distinct v2 JSON
schema. The content-addressed identity is stored once (default cap 16 MiB).
Every load/publish rehashes every original dependency, including resolving its
original alias anew. A manifest hash alone is never transitive verification.

v1 remains a separate explicit reader/writer. Mixed chains are refused. No
automatic migration, repair, optimizer acceptance, or budget reset is provided.
The private-directory, immutable-external-files, POSIX and caller-accounting
contracts of source_checkpoint_store apply unchanged.
"""
import hashlib
import os
from pathlib import Path

from .source_checkpoint_store import (
    StageCheckpointStore as _LegacyStore, StoredStage, _canonical, _fingerprints,
    _integer, _regular, _verify_artifacts, _FRAME, _MAGIC, _HEADER_LIMIT, _FIELDS)

_SCHEMA = 'source-stage-storage-v2'
_V2_FIELDS = (_FIELDS - {'identity'}) | {'identity_manifest'}


class StageCheckpointStore(_LegacyStore):
    def __init__(self, directory, identity_files, *, maximum_wall_ns,
                 maximum_trace_bytes, maximum_payload_bytes=16*1024*1024,
                 maximum_identity_bytes=16*1024*1024):
        # Reject oversized real provenance before allocating the store or running
        # any native work. This is a separate bound, not a larger frame header.
        limit = _integer(maximum_identity_bytes, 1, 64*1024*1024)
        identity = _fingerprints(identity_files)
        encoded = _canonical(dict(schema='source-stage-identity-v2', files=identity))
        if len(encoded) > limit:
            raise ValueError('IDENTITY_MANIFEST_LIMIT')
        super().__init__(directory, identity_files, maximum_wall_ns=maximum_wall_ns,
            maximum_trace_bytes=maximum_trace_bytes, maximum_payload_bytes=maximum_payload_bytes)
        if self._identity != identity:
            raise ValueError('STALE_CHECKPOINT_IDENTITY')
        self._limits['maximum_identity_bytes'] = limit
        self._manifest_bytes = encoded
        digest = hashlib.sha256(encoded).hexdigest()
        self._manifest_reference = dict(name=f'identity-{digest}.json', sha256=digest,
            size=len(encoded), file_count=len(identity))
        self.preflight()

    def _header(self, stage, previous, wall, trace, artifacts):
        return dict(schema=_SCHEMA, stage=stage, next_stage=stage+1,
            parent_sha256=previous.sha256 if previous else None,
            identity_manifest=self._manifest_reference, limits=self._limits,
            cumulative_wall_ns=wall, cumulative_trace_bytes=trace, artifacts=artifacts)

    def preflight(self, *, artifact_paths=None):
        """Check real identity and worst-width header before a caller starts work.

        Future artifact paths need not exist. Callers must supply their complete
        maximum artifact set (e.g. all planned pass traces) for an admission
        guarantee. Paths/aliases must remain fixed; sizes must fit signed 64-bit.
        With no paths this only checks the identity/base header, not future traces.
        Publish still verifies actual files, sizes and the actual header limit.
        """
        self._identity_current()
        paths = {} if artifact_paths is None else artifact_paths
        if not isinstance(paths, dict): raise ValueError('CHECKPOINT_FILE_MAP_REQUIRED')
        artifacts = {}
        for label, path in paths.items():
            if type(label) is not str or not label: raise ValueError('CHECKPOINT_FILE_LABEL_REQUIRED')
            artifacts[label] = dict(path=str(Path(path).resolve()), size=2**63-1, sha256='f'*64)
        h = self._header(99999999, None, self._limits['maximum_wall_ns'],
                         self._limits['maximum_trace_bytes'], artifacts)
        h['parent_sha256'] = 'f'*64
        size = len(_canonical(h))
        if size > _HEADER_LIMIT: raise ValueError('CHECKPOINT_HEADER_LIMIT')
        return dict(schema=_SCHEMA, identity_file_count=len(self._identity),
            identity_manifest_bytes=len(self._manifest_bytes),
            maximum_identity_bytes=self._limits['maximum_identity_bytes'],
            header_upper_bound_bytes=size, maximum_header_bytes=_HEADER_LIMIT,
            planned_artifact_count=len(artifacts), scientific_validation='NOT_GRANTED')

    def _read_manifest(self, reference):
        if _canonical(reference) != _canonical(self._manifest_reference):
            raise ValueError('STALE_CHECKPOINT_IDENTITY_MANIFEST')
        path = self.directory / self._manifest_reference['name']
        with _regular(path) as stream:
            if os.fstat(stream.fileno()).st_size != len(self._manifest_bytes):
                raise ValueError('IDENTITY_MANIFEST_SIZE')
            data = stream.read(len(self._manifest_bytes)+1)
        # Exact canonical bytes bind schema, labels, paths, sizes, hashes and count.
        # Original files were independently rehashed by _latest_locked/preflight.
        if data != self._manifest_bytes: raise ValueError('IDENTITY_MANIFEST_MISMATCH')

    def _ensure_manifest(self):
        """Exclusive lock held. Never replace or repair an existing object."""
        path = self.directory / self._manifest_reference['name']
        if not os.path.lexists(path):
            temporary = self._write_pending(self._manifest_bytes)
            self._identity_current()
            os.link(temporary, path)
        self._read_manifest(self._manifest_reference)
        # Sync even an orphan manifest from an interrupted earlier publication.
        fd = os.open(self.directory, os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
        try: os.fsync(fd)
        finally: os.close(fd)

    def _validate_header(self, header, index, previous):
        if (not isinstance(header, dict) or set(header) != _V2_FIELDS or
                header['schema'] != _SCHEMA):
            raise ValueError('INVALID_CHECKPOINT_HEADER')
        self._read_manifest(header['identity_manifest'])
        # Reuse all v1 sequence, parent, budget and artifact validation. Only the
        # on-disk identity representation differs; the full expected identity is
        # still constructed from real files, never trusted from the sidecar.
        translated = dict(header)
        del translated['identity_manifest']
        translated.update(schema='source-stage-storage-v1', identity=self._identity)
        super()._validate_header(translated, index, previous)

    def publish(self, stage, *, payload, cumulative_wall_ns, cumulative_trace_bytes, artifacts):
        _integer(stage, 0, 99999999)
        if type(payload) is not bytes or not 0 < len(payload) <= self._limits['maximum_payload_bytes']:
            raise ValueError('INVALID_CHECKPOINT_PAYLOAD')
        with self._lock(True):
            previous = self._latest_locked()
            if stage != (previous.next_stage if previous else 0):
                raise ValueError('CHECKPOINT_STAGE_SEQUENCE')
            records = _fingerprints(artifacts)
            for record in records.values(): _integer(record['size'], 0, 2**63-1)
            h = self._header(stage, previous, cumulative_wall_ns, cumulative_trace_bytes, records)
            header = _canonical(h)
            if len(header) > _HEADER_LIMIT: raise ValueError('CHECKPOINT_HEADER_LIMIT')
            self._ensure_manifest()
            self._validate_header(h, stage, previous)
            body = _FRAME.pack(_MAGIC, len(header), len(payload)) + header + payload
            digest = hashlib.sha256(body).hexdigest()
            target = self.directory / f'stage-{stage:08d}.ckpt'
            receipt = self.directory / f'commit-{stage:08d}.sha256'
            fd = os.open(self.directory, os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
            linked_stage = linked_receipt = False
            try:
                temporary = self._write_pending(body + digest.encode())
                self._identity_current()
                self._read_manifest(self._manifest_reference)
                _verify_artifacts(records)
                os.link(temporary, target); linked_stage = True
                os.fsync(fd)
                pending_receipt = self._write_pending(f'{stage:08d} {digest}\n'.encode())
                os.link(pending_receipt, receipt); linked_receipt = True
                os.fsync(fd)
            except BaseException:
                if linked_stage or linked_receipt:
                    if linked_receipt: os.unlink(receipt)
                    if linked_stage: os.unlink(target)
                    os.fsync(fd)
                raise
            finally: os.close(fd)
            return StoredStage(stage, stage+1, payload, digest, h['parent_sha256'],
                cumulative_wall_ns, cumulative_trace_bytes, target)
