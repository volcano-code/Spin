"""Explicit contact scenarios and record coverage; no geometry or physics grant.

The reconstructed object's spin rate is an input, not a paper parameter or a
universal stability threshold. An actual simulator must populate the records.
"""
from dataclasses import dataclass
import math
import re

PROFILES = ('default', 'robust')
STEP_DENOMINATORS = (240, 480, 960)


@dataclass(frozen=True)
class ContactScenario:
    repeat: int
    profile: str
    step_denominator: int
    direction: int | None

    @property
    def key(self):
        return (self.repeat, self.profile, self.step_denominator, self.direction)

    @property
    def tilt(self):
        return 0.0 if self.direction is None else 2.0

    @property
    def azimuth(self):
        return 0.0 if self.direction is None else self.direction * math.pi / 4

    @property
    def dt(self):
        return 1.0 / self.step_denominator


def scenarios():
    return tuple(ContactScenario(repeat, profile, step, direction)
                 for repeat in (0, 1) for profile in PROFILES
                 for step in STEP_DENOMINATORS for direction in (None, *range(8)))


def _finite_tree(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError('NONFINITE_RECORD')
    if isinstance(value, dict):
        for item in value.values():
            _finite_tree(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_tree(item)


def _near(actual, expected, label):
    if type(actual) not in (int, float) or not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
        raise ValueError('PROTOCOL_MISMATCH_' + label)


def audit_records(records, *, artifact_digest, spin_hz):
    """Reject inconsistent records, report missing scenarios, never grant stability.

    Directions are None for zero-tilt nominal and integers 0..7 for 2-degree
    perturbations. A direction-0 perturbation cannot stand in for nominal.
    """
    if not isinstance(artifact_digest, str) or not re.fullmatch('[0-9a-f]{64}', artifact_digest):
        raise ValueError('ARTIFACT_DIGEST_REQUIRED')
    if type(spin_hz) not in (int, float) or not math.isfinite(spin_hz) or spin_hz <= 0:
        raise ValueError('INVALID_SPIN_RATE')
    expected = {s.key: s for s in scenarios()}
    observed = {}
    for record in records:
        if not isinstance(record, dict) or 'direction' not in record:
            raise ValueError('EXPLICIT_SCENARIO_ID_REQUIRED')
        _finite_tree(record)
        repeat, profile = record.get('repeat'), record.get('profile')
        direction, denominator = record.get('direction'), record.get('step_denominator')
        if type(repeat) is not int or type(profile) is not str or type(denominator) is not int or (direction is not None and type(direction) is not int):
            raise ValueError('EXPLICIT_SCENARIO_ID_REQUIRED')
        key = (repeat, profile, denominator, direction)
        if key not in expected:
            raise ValueError('UNREGISTERED_SCENARIO')
        if key in observed:
            raise ValueError('DUPLICATE_SCENARIO')
        scenario = expected[key]
        if record.get('artifact_digest') != artifact_digest:
            raise ValueError('ARTIFACT_MISMATCH')
        if record.get('protocol') != 'PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT' or type(record.get('constraint_count')) is not int or record['constraint_count'] != 0:
            raise ValueError('FREE_CONTACT_PROTOCOL_REQUIRED')
        for field, value in (('initial_tilt', scenario.tilt), ('azimuth', scenario.azimuth),
                             ('dt', scenario.dt), ('duration', 8.0), ('hz', spin_hz),
                             ('friction', 0.8)):
            _near(record.get(field), value, field)
        samples = record.get('samples')
        if not isinstance(samples, list) or len(samples) < 2:
            raise ValueError('TRAJECTORY_REQUIRED')
        times = [s.get('time') for s in samples]
        if any(type(t) not in (float,int) or not math.isfinite(t) for t in times):
            raise ValueError('INVALID_SAMPLE_TIMES')
        _near(times[0], 0.0, 'sample_start')
        _near(times[-1], 8.0, 'sample_end')
        if any(b <= a for a,b in zip(times,times[1:])):
            raise ValueError('NONMONOTONIC_TRAJECTORY')
        if any(abs(t*denominator-round(t*denominator)) > 1e-8 for t in times):
            raise ValueError('SAMPLES_NOT_ON_RECORDED_TIMESTEP')
        _near(samples[0].get('tilt_degrees'), scenario.tilt, 'initial_observed_tilt')
        observed[key] = record
    missing = [s.key for s in scenarios() if s.key not in observed]
    return dict(scope='RECORD_IDENTITY_AND_PROTOCOL_COVERAGE_ONLY',
                status='INCOMPLETE' if missing else 'COMPLETE_PROTOCOL_COVERAGE',
                observed_scenarios=len(observed), expected_scenarios=len(expected),
                observed_nominal=sum(k[3] is None for k in observed), missing_scenarios=missing,
                geometry_acceptance='NOT_GRANTED', convergence='NOT_EVALUATED',
                determinism='NOT_EVALUATED', physics_acceptance='NOT_GRANTED',
                paper_reproduction='NOT_GRANTED')
