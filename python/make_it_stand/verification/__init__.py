"""Independent, standard-library geometry verification foundations.

The checker exposed here validates finite, consistently oriented closed
triangle 2-manifolds before integrating signed tetrahedra.  Reversed closed
components can represent cavities.  This is not a self-intersection,
containment, wall-thickness, or material-overlap certification.
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
]
