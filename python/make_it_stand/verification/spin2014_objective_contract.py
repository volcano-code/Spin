"""Immutable original-objective binding for the shared material pipeline.

This contract supplies independent exact equations, not optimization or a
geometry certificate. Production Newton must consume the native 2014 adapter.
"""
from dataclasses import asdict, dataclass
import hashlib
import json

from ..spin import Spin2014QuadricObjectiveConfig
from ._spin2014_authority import independent_authority_physics


@dataclass(frozen=True)
class Spin2014ObjectiveContract:
    settings: Spin2014QuadricObjectiveConfig

    def __post_init__(self):
        if type(self.settings) is not Spin2014QuadricObjectiveConfig:
            raise TypeError('expected immutable Spin2014QuadricObjectiveConfig')

    @property
    def mode(self):
        return {'YOYO': 'YO_YO', 'STATIC_LOWERING': 'STATIC_BALANCE'}.get(
            self.settings.objective, 'SPINNING_TOP')

    @property
    def definition(self):
        return dict(schema='spin2014-shared-objective-contract-v1',
                    objective_profile='SPIN2014_OBJECTIVES_QUADRIC_V1',
                    settings=asdict(self.settings), mode=self.mode,
                    density_convention='MOMENTS_ALREADY_DENSITY_WEIGHTED',
                    moment_order='M,Lx,Ly,Lz,Qxx,Qyy,Qzz,Qyz,Qxz,Qxy',
                    production_precision='BINARY64_COEFFICIENTS_EXACT_GEOMETRY_V1',
                    reference_precision='EXACT_RATIONAL_V1')

    @property
    def identity(self):
        encoded = json.dumps(self.definition, sort_keys=True, separators=(',', ':'),
                             allow_nan=False).encode()
        return hashlib.sha256(encoded).hexdigest()

    def require_mode(self, mode):
        if mode != self.mode:
            raise ValueError('original objective and request mode disagree')

    def independent_reference(self, raw_moments):
        """Evaluate canonical rational moments without multiplying density again."""
        return independent_authority_physics(raw_moments, asdict(self.settings))
