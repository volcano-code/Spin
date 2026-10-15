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
from .derived_material import derive_float_material
from .feature_thickness_contract import bind_feature_thickness_source_map
from .faster_material import (
    FasterMaterialResult,
    FloatingIslandChainResult,
    FrozenFasterMaterialComposition,
    compose_sealed_faster_material_as_fixed_region,
    export_spin_it_faster_quadric_material_step,
    quadric_coefficients_to_cad_frame,
    recompute_sealed_buoyant_stability_at_pose,
    seal_buoyant_plane_result_as_artifact,
    solve_floating_island_b_from_sealed_a,
    solve_floating_island_b_buoyant_from_sealed_a,
    solve_cherchi_geometry_and_seal_faster,
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
    "SpinArtifactProblem",
    "SealedSpinArtifact",
    "save_spin_artifact",
    "load_spin_artifact",
    "derive_float_material",
    "bind_feature_thickness_source_map",
    "FasterMaterialResult",
    "FloatingIslandChainResult",
    "FrozenFasterMaterialComposition",
    "compose_sealed_faster_material_as_fixed_region",
    "export_spin_it_faster_quadric_material_step",
    "quadric_coefficients_to_cad_frame",
    "recompute_sealed_buoyant_stability_at_pose",
    "seal_buoyant_plane_result_as_artifact",
    "solve_floating_island_b_from_sealed_a",
    "solve_floating_island_b_buoyant_from_sealed_a",
    "solve_cherchi_geometry_and_seal_faster",
]
