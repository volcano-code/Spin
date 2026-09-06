"""Independent geometry verification foundations.

The checker exposed here validates finite, consistently oriented closed
triangle 2-manifolds before integrating signed tetrahedra.  Reversed closed
components can represent cavities.  This is not a self-intersection,
containment, wall-thickness, or material-overlap certification.

The separate spin_problem/spin_artifact/spin_bundle modules add the strict
native exact HOLLOW_ONLY_MESH_MATERIAL_V1 path. Their native imports are lazy.
"""

from .mass_moments import (
    DensityRegion,
    GeometricMoments,
    GeometryValidationError,
    MassProperties,
    Matrix3,
    TriangleMesh,
    Vector3,
    aggregate_mass_properties,
    compute_mass_properties,
    integrate_geometric_moments,
)
from .spin_problem import SpinArtifactProblem
from .spin_artifact import SealedSpinArtifact
from .spin_bundle import save_spin_artifact, load_spin_artifact

__all__ = [
    "DensityRegion",
    "GeometricMoments",
    "GeometryValidationError",
    "MassProperties",
    "Matrix3",
    "TriangleMesh",
    "Vector3",
    "aggregate_mass_properties",
    "compute_mass_properties",
    "integrate_geometric_moments",
    "SpinArtifactProblem",
    "SealedSpinArtifact",
    "save_spin_artifact",
    "load_spin_artifact",
]
