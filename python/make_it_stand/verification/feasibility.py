"""Independent exact certificates for a represented discrete moment LP.

Each finite binary floating input is interpreted as its exact integer ratio.
For a CCW edge (dx, dy), a=(dy,-dx), b=a.vertex, C_j=a.L_j-b*m_j,
d=b*m_fixed-a.L_fixed. All subsequent arithmetic uses Fraction. This is
not the native LP's rounded unit-normal matrix, nor continuous geometric
carving: x in [0,1] retains fractions on the specified discrete fixed cells.

For nonnegative nonzero lambda, every boxed x has
lambda.C.x >= sum_j min(0, (C.T.lambda)_j). A strict excess of this lower
bound over lambda.d proves infeasibility of this represented system only.
No solver status, tolerance, normalization, or stationarity test is proof.
Input hashes bind numeric inputs, ordered IDs and target; they do not prove
where those numbers came from or attest geometry, source, or a loaded build.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from fractions import Fraction as Q
import hashlib
import importlib
import json
import math
from pathlib import Path
from typing import Any

STATUSES = frozenset({'FEASIBLE_WITNESS', 'INFEASIBLE_CERTIFIED', 'UNKNOWN'})
SCOPE = 'REPRESENTED_DISCRETE_MOMENT_LP'


def _number(value: Any, *, candidate: bool = False) -> Q:
    if type(value) is int:
        return Q(value)
    if type(value) is float and math.isfinite(value):
        return Q(*value.as_integer_ratio())
    if candidate and type(value) is Q:
        return value
    raise ValueError('expected finite binary number (not bool/string)')


def _vector(values: Any, size: int, *, candidate: bool = False) -> tuple[Q, ...]:
    if isinstance(values, (str, bytes, dict)):
        raise ValueError('expected numeric sequence')
    result = tuple(_number(v, candidate=candidate) for v in values)
    if len(result) != size:
        raise ValueError(f'expected {size} entries, received {len(result)}')
    return result


def _encode(value: Any) -> Any:
    if isinstance(value, Q):
        return [value.numerator, value.denominator]
    if isinstance(value, dict):
        return {k: _encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(v) for v in value]
    return value


def _hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(_encode(value), sort_keys=True,
                                    separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class MomentCell:
    cell_id: int
    mass: float
    first_moment: tuple[float, float]


@dataclass(frozen=True)
class MomentSnapshot:
    fixed_mass: float
    fixed_first_moment: tuple[float, float]
    cells: tuple[MomentCell, ...]
    polygon: tuple[tuple[float, float], ...]
    tolerance: float

    @classmethod
    def from_native(cls, value: Any) -> MomentSnapshot:
        # Native numpy arrays are copied to builtin finite scalars. The strict
        # entry point checks native types before reaching this adapter.
        return cls(value.fixed_mass, tuple(value.fixed_first_moment_2d.tolist()),
                   tuple(MomentCell(c.cell_id, c.mass, tuple(c.first_moment_2d.tolist()))
                         for c in value.variable_cells),
                   tuple(tuple(row) for row in value.safe_polygon_2d.tolist()),
                   value.feasibility_tolerance)


@dataclass(frozen=True)
class ExactSystem:
    C: tuple[tuple[Q, ...], ...]
    d: tuple[Q, ...]
    fixed_mass: Q
    fixed_first_moment: tuple[Q, Q]
    masses: tuple[Q, ...]
    moments: tuple[tuple[Q, Q], ...]
    input_hash: str
    system_hash: str


def build_system(snapshot: MomentSnapshot) -> ExactSystem:
    """Validate and rebuild unnormalized moment rows without native helpers.

    Every other vertex must lie on/left of every oriented edge; positive
    signed area, distinct vertices and forward collinear edges ensure a
    simple convex CCW boundary. Redundant collinear vertices are accepted.
    """
    if not isinstance(snapshot, MomentSnapshot):
        raise ValueError('expected MomentSnapshot, not historical statistics')
    mass = _number(snapshot.fixed_mass)
    fixed = _vector(snapshot.fixed_first_moment, 2)
    tolerance = _number(snapshot.tolerance)
    if mass <= 0 or tolerance <= 0:
        raise ValueError('fixed mass and tolerance must be positive')
    polygon = tuple(_vector(v, 2) for v in snapshot.polygon)
    if len(polygon) < 3 or len(set(polygon)) != len(polygon):
        raise ValueError('polygon requires >=3 distinct vertices, without repeats')
    edges = tuple((polygon[(i+1) % len(polygon)][0]-p[0],
                   polygon[(i+1) % len(polygon)][1]-p[1]) for i, p in enumerate(polygon))
    area = sum(p[0]*polygon[(i+1) % len(polygon)][1] -
               p[1]*polygon[(i+1) % len(polygon)][0] for i, p in enumerate(polygon))
    if area <= 0:
        raise ValueError('polygon requires positive CCW area')
    for i, ((dx, dy), p) in enumerate(zip(edges, polygon)):
        if dx == dy == 0:
            raise ValueError('zero-length polygon edge')
        if any(dx*(q[1]-p[1])-dy*(q[0]-p[0]) < 0 for q in polygon):
            raise ValueError('polygon is not a simple convex CCW boundary')
        ex, ey = edges[(i+1) % len(edges)]
        if dx*ey-dy*ex == 0 and dx*ex+dy*ey <= 0:
            raise ValueError('overlapping/backtracking collinear edges')
    cells = []
    seen = set()
    for cell in snapshot.cells:
        if not isinstance(cell, MomentCell) or type(cell.cell_id) is not int or cell.cell_id < 0 or cell.cell_id in seen:
            raise ValueError('cell IDs must be unique nonnegative integers')
        seen.add(cell.cell_id)
        m = _number(cell.mass)
        moment = _vector(cell.first_moment, 2)
        if m < 0 or (m == 0 and any(moment)):
            raise ValueError('negative mass or zero-mass cell with nonzero moments')
        cells.append((cell.cell_id, m, moment))
    masses = tuple(c[1] for c in cells)
    moments = tuple(c[2] for c in cells)
    if mass + sum(masses) <= 0:
        raise ValueError('positive total mass required')
    C, d = [], []
    for (dx, dy), p in zip(edges, polygon):
        ax, ay = dy, -dx
        b = ax*p[0]+ay*p[1]
        C.append(tuple(ax*l[0]+ay*l[1]-b*m for m, l in zip(masses, moments)))
        d.append(b*mass-ax*fixed[0]-ay*fixed[1])
    identity = _hash(('represented-binary-moment-input-v1', mass, fixed, cells, polygon, tolerance))
    system_hash = _hash(('exact-bounded-moment-system-v1', identity, C, d, (0, 1)))
    return ExactSystem(tuple(C), tuple(d), mass, fixed, masses, moments, identity, system_hash)


@dataclass(frozen=True)
class ExactCertificate:
    reason: str = 'MISSING_CANDIDATE'
    weights: tuple[Q, ...] = ()
    column_sums: tuple[Q, ...] = ()
    lower_bound: Q | None = None
    weighted_rhs: Q | None = None
    margin: Q | None = None

    @property
    def valid(self) -> bool:
        return self.reason == 'STRICT_POSITIVE_MARGIN' and self.margin is not None and self.margin > 0


def _certificate(system: ExactSystem, weights: Any) -> ExactCertificate:
    if weights is None:
        return ExactCertificate()
    try:
        weights = _vector(weights, len(system.d), candidate=True)
        if any(w < 0 for w in weights) or not any(weights):
            raise ValueError('weights must be nonnegative and nonzero')
    except (TypeError, ValueError, OverflowError) as error:
        return ExactCertificate('INVALID_CANDIDATE: ' + str(error))
    sums = tuple(sum((row[j]*w for row, w in zip(system.C, weights)), Q(0)) for j in range(len(system.masses)))
    lower = sum((min(Q(0), s) for s in sums), Q(0))
    rhs = sum((w*d for w, d in zip(weights, system.d)), Q(0))
    margin = lower-rhs
    return ExactCertificate('STRICT_POSITIVE_MARGIN' if margin > 0 else 'NONPOSITIVE_EXACT_MARGIN',
                            weights, sums, lower, rhs, margin)


@dataclass(frozen=True)
class AlgebraicWitness:
    reason: str = 'MISSING_CANDIDATE'
    retained_fractions: tuple[Q, ...] = ()
    row_slacks: tuple[Q, ...] = ()
    mass: Q | None = None
    first_moment: tuple[Q, ...] = ()
    center: tuple[Q, ...] = ()

    @property
    def valid(self) -> bool:
        return self.reason == 'EXACT_LP_WITNESS_ONLY'


def _witness(system: ExactSystem, retained: Any) -> AlgebraicWitness:
    if retained is None:
        return AlgebraicWitness()
    try:
        x = _vector(retained, len(system.masses), candidate=True)
        if any(v < 0 or v > 1 for v in x):
            raise ValueError('retained fractions outside [0,1]')
    except (TypeError, ValueError, OverflowError) as error:
        return AlgebraicWitness('INVALID_CANDIDATE: ' + str(error))
    slacks = tuple(d-sum((c*v for c, v in zip(row, x)), Q(0)) for row, d in zip(system.C, system.d))
    mass = system.fixed_mass + sum((m*v for m, v in zip(system.masses, x)), Q(0))
    moments = tuple(system.fixed_first_moment[k] + sum((l[k]*v for l, v in zip(system.moments, x)), Q(0)) for k in range(2))
    return AlgebraicWitness('EXACT_LP_WITNESS_ONLY' if all(s >= 0 for s in slacks) else 'EXACT_ROW_VIOLATION',
                            x, slacks, mass, moments, tuple(l/mass for l in moments))


@dataclass(frozen=True)
class ObservedNativeIdentity:
    path: str = ''
    sha256: str = ''
    observation: str = 'NOT_OBSERVED'


@dataclass(frozen=True)
class FeasibilityReport:
    reason: str
    input_hash: str = ''
    system_hash: str = ''
    certificate: ExactCertificate = field(default_factory=ExactCertificate)
    witness: AlgebraicWitness = field(default_factory=AlgebraicWitness)
    checker_sha256: str = ''
    native_identity: ObservedNativeIdentity = field(default_factory=ObservedNativeIdentity)
    native_input_hash: str = ''
    native_pre_input_hash: str = ''
    native_post_input_hash: str = ''
    native_status: str = ''
    native_primal_status: str = ''
    native_certificate_status: str = ''
    candidate_origin: str = 'SUPPLIED_EXACT_ROW_WEIGHTS'
    scope: str = field(default=SCOPE, init=False)
    geometry_validation: str = field(default='NOT_VERIFIED', init=False)
    continuous_geometry: str = field(default='UNKNOWN', init=False)
    physical_validation: str = field(default='NOT_IN_SCOPE', init=False)
    build_source_attestation: str = field(default='NOT_VERIFIED', init=False)

    @property
    def status(self) -> str:
        # No public geometry-valid flag or path to FEASIBLE_WITNESS. That status
        # requires future SealedSpinArtifact auditing outside this checker.
        return 'INFEASIBLE_CERTIFIED' if self.reason == 'EXACT_INFEASIBILITY_CERTIFICATE' and self.certificate.valid else 'UNKNOWN'

    def to_dict(self) -> dict[str, Any]:
        result = _encode(asdict(self))
        result.update(schema='represented-discrete-feasibility-v1', status=self.status)
        result['certificate']['valid'] = self.certificate.valid
        result['witness']['valid'] = self.witness.valid
        return result


def verify_snapshot(snapshot: MomentSnapshot, *, weights: Any = None,
                    retained_fractions: Any = None, expected_input_hash: str | None = None) -> FeasibilityReport:
    """Check algebra on supplied moments; this entry does not attest native origin.

    expected_input_hash is only a consistency guard, never proof of construction.
    Encoded rational quantities are [numerator, denominator] integer pairs.
    """
    checker = _file_hash(Path(__file__))
    try:
        system = build_system(snapshot)
    except (TypeError, ValueError, AttributeError, OverflowError) as error:
        return FeasibilityReport('INVALID_INPUT: ' + str(error), checker_sha256=checker)
    report = FeasibilityReport('', system.input_hash, system.system_hash, checker_sha256=checker)
    if expected_input_hash is not None and expected_input_hash != system.input_hash:
        return replace(report, reason='INPUT_IDENTITY_MISMATCH')
    certificate = _certificate(system, weights)
    witness = _witness(system, retained_fractions)
    reason = 'EXACT_INFEASIBILITY_CERTIFICATE' if certificate.valid else certificate.reason
    return replace(report, reason=reason, certificate=certificate, witness=witness)


def _native_statistics(result: Any, rows: int, columns: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Reject malformed observations; the booleans/residuals cannot certify."""
    for count in (result.input_cell_count, result.variable_cell_count, result.certificate_system_nonzero_count):
        if type(count) is not int or count < 0:
            raise ValueError('invalid native count')
    if result.input_cell_count != columns or result.variable_cell_count != columns:
        raise ValueError('native cell count differs from snapshot')
    if result.status not in ('ALREADY_SAFE', 'RELAXATION_REACHABLE', 'PROVEN_UNREACHABLE', 'INDETERMINATE'):
        raise ValueError('invalid native status')
    if not isinstance(result.primal_solver_status, str) or not result.primal_solver_status:
        raise ValueError('empty/malformed native primal status')
    if not isinstance(result.certificate_solver_status, str) or type(result.certificate.valid) is not bool:
        raise ValueError('malformed native certificate statistics')
    for value in (result.fixed_mass, result.final_mass, result.maximum_primal_violation,
                  result.certificate.minimum_weight, result.certificate.normalization_residual,
                  result.certificate.stationarity_infinity_norm, result.certificate.contradiction):
        _number(value)
    _vector(result.final_first_moment_2d.tolist(), 2)
    _vector(result.final_center_2d.tolist(), 2)
    weights = tuple(result.certificate.weights.tolist())
    retained = tuple(result.retained_fractions.tolist())
    if weights:
        if not result.certificate_solver_status:
            raise ValueError('empty native certificate solver status with candidate')
        _vector(weights, rows+2*columns)
        if any(w < 0 for w in weights) or not any(weights):
            raise ValueError('native candidate has negative or zero weights')
    if retained:
        _vector(retained, columns)
    return weights, retained


def check_native_feasibility(subject: Any, standing_index: int = 0) -> FeasibilityReport:
    """Observe pre-input, analyze, observe post-input, independently verify.

    Accepts a real native BalanceModel or immutable typed fixture input. This
    sequence detects differing numeric identities, but does not promise
    atomicity, detect ABA geometry changes, or attest in-memory code. File
    hashes observe disk bytes only; build/source attestation stays unresolved.
    Native polygon weights are used unchanged on unnormalized exact rows:
    they are fresh candidates, not an assertion of equality to the native
    unit-normal proof. No rescaling or clamping is performed.
    """
    report = FeasibilityReport('NATIVE_CHECK_NOT_COMPLETED', checker_sha256=_file_hash(Path(__file__)),
                               candidate_origin='NATIVE_POLYGON_PREFIX_UNCHANGED_ON_EXACT_UNNORMALIZED_ROWS')
    try:
        native = importlib.import_module('_make_it_stand')
        if not all(hasattr(native, name) for name in ('CarvingReachabilityInput', 'BalanceModel', 'analyze_carving_reachability_input')):
            return replace(report, reason='MISSING_NATIVE_REACHABILITY_BACKEND')
        path = Path(native.__file__).resolve(strict=True)
        report = replace(report, native_identity=ObservedNativeIdentity(str(path), _file_hash(path), 'OBSERVED_DISK_BYTES_NOT_IN_MEMORY_ATTESTATION'))
    except (ImportError, OSError, TypeError, AttributeError) as error:
        return replace(report, reason='MISSING_NATIVE_BACKEND: ' + str(error))
    if type(subject) not in (native.BalanceModel, native.CarvingReachabilityInput):
        return replace(report, reason='INVALID_NATIVE_SUBJECT')
    if type(standing_index) is not int or standing_index < 0:
        return replace(report, reason='INVALID_STANDING_INDEX')
    is_model = type(subject) is native.BalanceModel
    try:
        before = subject.carving_reachability_snapshot(standing_index) if is_model else subject
        data = before.input if is_model else before
        snapshot = MomentSnapshot.from_native(data)
        report = replace(report, native_pre_input_hash=before.input_hash)
        # Independent input validation happens even for native early returns.
        algebra = verify_snapshot(snapshot)
        report = replace(report, input_hash=algebra.input_hash, system_hash=algebra.system_hash)
        result = subject.analyze_carving_reachability(standing_index) if is_model else native.analyze_carving_reachability_input(subject)
        if type(result) is not native.CarvingReachabilityResult:
            return replace(report, reason='INVALID_NATIVE_RESULT_TYPE')
        report = replace(report, native_input_hash=result.input_hash, native_status=result.status,
                         native_primal_status=result.primal_solver_status,
                         native_certificate_status=result.certificate_solver_status)
        after = subject.carving_reachability_snapshot(standing_index) if is_model else subject
        post = MomentSnapshot.from_native(after.input if is_model else after)
        report = replace(report, native_post_input_hash=after.input_hash)
        if (before.input_hash != result.input_hash or result.input_hash != after.input_hash or snapshot != post
                or (is_model and before.preparation_status != after.preparation_status)):
            return replace(report, reason='INPUT_IDENTITY_MISMATCH')
        if is_model and before.preparation_status != 'READY':
            return replace(report, reason=before.preparation_status + '; ' + algebra.reason)
        if algebra.reason.startswith('INVALID_INPUT'):
            return replace(report, reason=algebra.reason)
        weights, retained = _native_statistics(result, len(snapshot.polygon), len(snapshot.cells))
        checked = verify_snapshot(snapshot, weights=weights[:len(snapshot.polygon)] if weights else None,
                                  retained_fractions=retained if retained or not snapshot.cells else None,
                                  expected_input_hash=report.input_hash)
        return replace(report, reason=checked.reason, certificate=checked.certificate, witness=checked.witness)
    except Exception as error:
        # Backend errors, failed snapshots and malformed observations are
        # explicit failures. Preserve any identities/statuses already observed.
        return replace(report, reason='NATIVE_CHECK_FAILURE: ' + type(error).__name__ + ': ' + str(error))
