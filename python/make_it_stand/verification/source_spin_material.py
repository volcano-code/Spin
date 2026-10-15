"""Actual source Spin2014 solve bound to independently reintegrated material.

Geometry stays in original model coordinates. Only the native moment claim is
converted to INPUT_NORMALIZED_V1. A saved spin_bundle authenticates geometry,
not solver execution; there is deliberately no trust-the-run-JSON loader.
All numerical and resource defaults are IMPLEMENTATION_DEFAULT.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import math
import numpy as np

from .spin_artifact import SealedSpinArtifact, _recheck_problem
from .mass_moments import TriangleMesh, compute_mass_properties
from .spin_common import (FRAME, ORDER, SCOPE, canonical, digest, native, parse, require,
                          scalar, vector, implementation_observation)
from .spin_problem import MeshSnapshot, SpinArtifactProblem


def _count(value, name, minimum=1, maximum=2**31-1):
    require(isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)),
            name + ': expected integer, bool is invalid')
    require(minimum <= int(value) <= maximum, name + ': outside supported range')
    return int(value)


def _numbers(obj, names, *, positive=False):
    for name in names.split():
        object.__setattr__(obj, name, scalar(getattr(obj, name), name, positive=positive))


@dataclass(frozen=True)
class SourceGeometrySettings:
    root_origin: tuple = (0., 0., 0.)
    root_side: float = 1.
    maximum_depth: int = 9
    maximum_nodes: int = 5000000

    def __post_init__(self):
        object.__setattr__(self, 'root_origin', tuple(vector(self.root_origin, 'root_origin')))
        _numbers(self, 'root_side', positive=True)
        object.__setattr__(self, 'maximum_depth', _count(self.maximum_depth, 'maximum_depth', 0, 20))
        object.__setattr__(self, 'maximum_nodes', _count(self.maximum_nodes, 'maximum_nodes', maximum=2**53))


@dataclass(frozen=True)
class AuditedSlqpSettings:
    maximum_iterations: int = 120
    maximum_line_search_steps: int = 24
    lbfgs_memory: int = 8
    equality_tolerance: float = 1e-8
    projected_gradient_tolerance: float = 1e-7
    active_set_tolerance: float = 1e-10
    merit_penalty: float = 100.
    line_search_reduction: float = .5
    armijo: float = 1e-4
    algorithm: str = 'AUDITED_KKT_RESTORATION_V1'
    restoration_trust_bound: float = 1.
    lp_iteration_limit: int = 10000
    lp_time_limit: float = 10.
    lp_feasibility_tolerance: float = 1e-9

    def __post_init__(self):
        require(self.algorithm == 'AUDITED_KKT_RESTORATION_V1', 'unsupported SLQP algorithm')
        for name in ('maximum_iterations','maximum_line_search_steps','lbfgs_memory','lp_iteration_limit'):
            object.__setattr__(self, name, _count(getattr(self, name), name))
        _numbers(self, 'equality_tolerance projected_gradient_tolerance merit_penalty line_search_reduction armijo restoration_trust_bound lp_time_limit lp_feasibility_tolerance', positive=True)
        _numbers(self, 'active_set_tolerance')
        require(0 <= self.active_set_tolerance < .5, 'invalid active_set_tolerance')
        require(self.line_search_reduction < 1 and self.armijo < 1, 'invalid line search control')
        require(1e-10 <= self.lp_feasibility_tolerance <= 1e-4, 'invalid LP feasibility tolerance')


_MODES = {'STATIC_LOWERING':'STATIC_BALANCE', 'YOYO':'YO_YO',
          'TOP_INERTIA_ONLY':'SPINNING_TOP', 'TOP_LOWERING_ONLY':'SPINNING_TOP',
          'TOP_MASS_AND_LOWERING':'SPINNING_TOP'}


@dataclass(frozen=True)
class SourceSolverSettings:
    objective: str = 'STATIC_LOWERING'
    formulation: str = 'PAPER_GIVENS'
    initial_level: int = 2
    maximum_optimization_level: int = 5
    maximum_adaptive_passes: int = 16
    maximum_leaf_count: int = 1000000
    maximum_output_cuboids: int = 1000000
    initial_fill: float = .5
    initial_phi: float = 0.
    fractional_epsilon: float = .05
    laplacian_weight: float = 1e-3
    inertia_weight: float = 1.
    center_weight: float = 1.
    binary_constraint_tolerance: float = 1e-6
    target_origin: tuple[float, float, float] = (0., 0., 0.)
    target_axes: tuple[tuple[float, float, float], ...] = (
        (1., 0., 0.), (0., 1., 0.), (0., 0., 1.),
    )
    slqp: AuditedSlqpSettings = field(default_factory=AuditedSlqpSettings)
    termination_policy: str = 'UNCHANGED_PARTITION_KKT_V1'

    def __post_init__(self):
        require(isinstance(self.objective, str) and self.objective in _MODES, 'invalid objective')
        require(self.formulation == 'PAPER_GIVENS', 'unsupported formulation')
        require(isinstance(self.termination_policy,str) and self.termination_policy in
                ('UNCHANGED_PARTITION_KKT_V1','PAPER_POST_TRANSITION_V1'),
                'unsupported termination policy')
        require(type(self.slqp) is AuditedSlqpSettings, 'expected immutable AuditedSlqpSettings')
        for name in ('initial_level', 'maximum_optimization_level'):
            object.__setattr__(self, name, _count(getattr(self, name), name, 0, 20))
        require(self.initial_level <= self.maximum_optimization_level, 'invalid optimization levels')
        for name in ('maximum_adaptive_passes','maximum_leaf_count','maximum_output_cuboids'):
            bound = (2**31-1)//6 if name == 'maximum_leaf_count' else 2**31-1
            object.__setattr__(self, name, _count(getattr(self, name), name, maximum=bound))
        _numbers(self, 'initial_fill initial_phi fractional_epsilon laplacian_weight inertia_weight center_weight')
        _numbers(self, 'binary_constraint_tolerance', positive=True)
        require(0 <= self.initial_fill <= 1 and 0 < self.fractional_epsilon < .5, 'invalid fill controls')
        require(min(self.laplacian_weight,self.inertia_weight,self.center_weight) >= 0, 'negative weight')
        require(self.objective != 'STATIC_LOWERING' or self.initial_phi == 0, 'static requires phi zero')
        origin=np.asarray(self.target_origin,dtype=float);axes=np.asarray(self.target_axes,dtype=float)
        require(origin.shape==(3,) and np.all(np.isfinite(origin)), 'invalid target origin')
        require(axes.shape==(3,3) and np.all(np.isfinite(axes)), 'invalid target axes')
        require(np.allclose(axes.T@axes,np.eye(3),rtol=0.,atol=1e-12) and
                math.isclose(float(np.linalg.det(axes)),1.,rel_tol=0.,abs_tol=1e-12),
                'target axes must be right-handed orthonormal')
        object.__setattr__(self,'target_origin',tuple(float(v) for v in origin))
        object.__setattr__(self,'target_axes',tuple(tuple(float(v) for v in row) for row in axes))


@dataclass(frozen=True)
class SourceBoundarySettings:
    maximum_input_cuboids: int = 10000
    maximum_output_vertices: int = 1000000
    maximum_output_triangles: int = 2000000
    maximum_intermediate_complexity: int = 4000000

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, _count(getattr(self, name), name, maximum=2**53))


@dataclass(frozen=True, init=False)
class SourceSpinMaterialResult:
    status: str
    artifact: SealedSpinArtifact | None
    seal_error: str | None
    scientific_digest: str
    _identity: str
    _outer: MeshSnapshot | None
    _material: MeshSnapshot | None
    _cavity: MeshSnapshot | None
    _boundary_material: MeshSnapshot | None
    _interior_material: MeshSnapshot | None
    _evidence: str
    _checkpoint: bytes | None

    def __init__(self):
        raise TypeError('use solve_and_seal_source_spin2014')

    @property
    def identity(self): return parse(self._identity)

    @property
    def evidence(self): return parse(self._evidence)

    @property
    def checkpoint(self): return self._checkpoint

    @property
    def outer(self): return None if self._outer is None else self._outer.arrays()

    @property
    def material(self): return None if self._material is None else self._material.arrays()

    @property
    def cavity(self): return None if self._cavity is None else self._cavity.arrays()

    @property
    def boundary_material(self):
        return None if self._boundary_material is None else self._boundary_material.arrays()

    @property
    def interior_material(self):
        return None if self._interior_material is None else self._interior_material.arrays()


def _normalize_model_moments(moments, normalization):
    """Full affine mass/first/second transformation with density/volume scaling."""
    m, x, y, z, xx, yy, zz, yz, xz, xy = moments
    o = np.array(normalization['origin'])
    b = np.array([x,y,z])
    q = np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
    q = q - np.outer(o,b) - np.outer(b,o) + m*np.outer(o,o)
    length, density = normalization['L'], normalization['R']
    first = (b-m*o)/(density*length**4)
    q /= density*length**5
    return [m/(density*length**3), *first.tolist(), *np.diag(q).tolist(),
            float(q[1,2]),float(q[0,2]),float(q[0,1])]


def _append_reversed(whole, cavity):
    if cavity.explicit_empty: return whole
    v,f = whole.arrays(); cv,cf = cavity.arrays()
    return MeshSnapshot((np.vstack((v,cv)),np.vstack((f,cf[:,::-1]+len(v)))))


def _independent_material_moments(mesh, density):
    """Independent signed-tetrahedron integral in source model coordinates."""
    require(not mesh.explicit_empty, 'cannot integrate an empty material component')
    vertices, faces = mesh.arrays()
    props = compute_mass_properties(TriangleMesh(
        tuple(map(tuple, vertices.tolist())), tuple(map(tuple, faces.tolist()))), density)
    q = props.second_moment
    return (props.mass, *props.first_moment, q[0][0], q[1][1], q[2][2],
            q[1][2], q[0][2], q[0][1])


def _dual_material_independent_audit(boundary_mesh, interior_mesh, boundary_density,
                                     interior_density, native_geometry, root_side):
    """Reintegrate each exported phase and gate sealing on agreement and partition."""
    boundary = np.asarray(_independent_material_moments(boundary_mesh, boundary_density))
    interior = np.asarray(_independent_material_moments(interior_mesh, interior_density))
    combined = boundary + interior
    claimed_boundary = np.asarray(native_geometry['boundary_material_moments'], dtype=float)
    claimed_interior = np.asarray(native_geometry['interior_material_moments'], dtype=float)
    claimed_combined = np.asarray(native_geometry['integrated_moments'], dtype=float)
    subtraction = np.asarray(native_geometry['source_subtraction_moments'], dtype=float)
    powers = np.asarray([3,4,4,4,5,5,5,5,5,5], dtype=float)
    scales = max(boundary_density, interior_density) * np.power(root_side, powers)
    errors = {
        'boundary_vs_native': float(np.max(np.abs(boundary-claimed_boundary)/scales)),
        'interior_vs_native': float(np.max(np.abs(interior-claimed_interior)/scales)),
        'sum_vs_native': float(np.max(np.abs(combined-claimed_combined)/scales)),
        'sum_vs_source_subtraction': float(np.max(np.abs(combined-subtraction)/scales)),
    }
    require(all(np.isfinite(value) and value <= 1e-10 for value in errors.values()),
            'dual-material independent mass-moment audit failed: '+canonical(errors))
    partition = native_geometry['dual_material_partition']
    require(partition['completed'] and partition['complete_partition'] and
            partition['no_positive_volume_overlap'],
            'dual-material exact partition audit failed: '+canonical(partition))
    return {'method':'INDEPENDENT_SIGNED_TETRAHEDRA_V1',
            'order':list(ORDER), 'density_applied_once':True,
            'boundary_moments':boundary.tolist(), 'interior_moments':interior.tolist(),
            'combined_moments':combined.tolist(), 'normalized_residuals':errors,
            'normalized_tolerance':1e-10, 'exact_partition':partition}


def _observation():
    from .run_evidence import _record
    result = implementation_observation()
    result['bridge'] = _record(Path(__file__))
    root = Path(__file__).resolve().parents[3]
    # Bounded file observations supplement the generic checker observations.
    # The runner captures the whole declared effective source policy separately.
    paths = ['bindings/SourceSpinMaterialBindings.cpp','bindings/SourceSpinMaterialBindings.h',
             'bindings/pybind_module.cpp','bindings/ExactMeshSequenceIngestion.h',
             'src/spin2014/SourceAdaptiveSpinSolver.cpp','src/spin2014/SourceMaterialBoundary.cpp',
             'src/spin2014/SourceSurfaceOctree.cpp','src/spin2014/SourceBinaryMoments.cpp']
    result['local_source'] = {name:_record(root/name) if (root/name).is_file() else None for name in paths}
    result['source_coverage'] = 'DECLARED_LOCAL_FILES_ONLY; null means unavailable'
    return result


def solve_and_seal_source_spin2014(problem, *, source_settings, solver_settings, boundary_settings,
                                  maximum_new_passes=None, resume_checkpoint=None,
                                  iteration_store_directory=None, iteration_store_max_bytes=1073741824):
    """One actual solve, original-coordinate export, and independent exact seal."""
    execution={}
    if maximum_new_passes is not None:
        execution['maximum_new_passes']=_count(maximum_new_passes,'maximum_new_passes')
        require(isinstance(iteration_store_directory,str) and Path(iteration_store_directory).is_absolute(),
                'staged material requires absolute iteration_store_directory')
    require(resume_checkpoint is None or (maximum_new_passes is not None and
            type(resume_checkpoint) is bytes and 0<len(resume_checkpoint)<=64*1024**2),
            'invalid resume checkpoint/control')
    if resume_checkpoint is not None:execution['resume_checkpoint']=resume_checkpoint
    if iteration_store_directory is not None:
        require(isinstance(iteration_store_directory,str) and '\0' not in iteration_store_directory,
                'invalid iteration_store_directory')
        execution['iteration_store_directory']=iteration_store_directory
        execution['iteration_store_max_bytes']=_count(iteration_store_max_bytes,'iteration_store_max_bytes',maximum=2**63-1)
    for value, cls in ((source_settings,SourceGeometrySettings),(solver_settings,SourceSolverSettings),
                       (boundary_settings,SourceBoundarySettings)):
        require(type(value) is cls, 'expected immutable '+cls.__name__)
        # Revalidate even if object.__setattr__ was used to bypass frozen setters.
        data = asdict(value)
        if cls is SourceSolverSettings: data['slqp'] = AuditedSlqpSettings(**data['slqp'])
        require(asdict(cls(**data)) == asdict(value), 'settings identity mismatch')
    require(solver_settings.maximum_optimization_level <= source_settings.maximum_depth,
            'optimization depth exceeds source geometry depth')
    require(isinstance(problem, SpinArtifactProblem), 'expected locked SpinArtifactProblem')
    problem = _recheck_problem(problem)
    require(problem.validation_scope == SCOPE and len(problem.fixed_regions) == 1 and len(problem.design_regions) == 1,
            'NOT_SUPPORTED: exactly one fixed shell and one design region in hollow-only scope required')
    boundary_density=problem.fixed_regions[0].density
    interior_density=problem.design_regions[0].density
    target = problem.target
    axes=np.asarray(solver_settings.target_axes,dtype=float)
    require(np.allclose(target['origin'],solver_settings.target_origin,rtol=0.,atol=1e-12) and
            np.allclose(target['axis'],axes[:,2],rtol=0.,atol=1e-12),
            'NOT_SUPPORTED: problem target must match source solver target frame')
    require(target['mode'] == _MODES[solver_settings.objective], 'NOT_SUPPORTED: objective/target mode mismatch')
    require(target['support_points'] == ([[0,0,0]] if target['mode']=='STATIC_BALANCE' else []),
            'NOT_SUPPORTED: unsupported support context')
    n = native()
    require(hasattr(n, 'solve_source_spin2014_material'), 'source bridge requires native 2014 plus exact geometry')
    before = _observation()
    solver_config=asdict(solver_settings)
    if boundary_density != interior_density:
        solver_config['boundary_density']=boundary_density
        solver_config['interior_density']=interior_density
    run = n.solve_source_spin2014_material(*problem.original.arrays(), *problem.design_regions[0].mesh.arrays(),
        source_settings=asdict(source_settings), solver_settings=solver_config,
        boundary_settings=asdict(boundary_settings), density=boundary_density, minimum_outer_wall=problem.minimum_outer_wall,**execution)
    diagnostics = run.diagnostics
    artifact, error = None, None
    outer = material = cavity = boundary_material = interior_material = None
    dual_material_audit = None
    if run.has_geometry:
        try:
            outer = MeshSnapshot(run.outer)
            require(outer == problem.original, 'native outer bytes differ from locked original')
            material = MeshSnapshot(run.material)
            cavity = MeshSnapshot(run.cavity if len(run.cavity[0]) else None)
            if boundary_density != interior_density:
                require(run.boundary_material is not None and run.interior_material is not None,
                        'native dual-material component boundaries are unavailable')
                boundary_material=MeshSnapshot(run.boundary_material if len(run.boundary_material[0]) else None)
                interior_material=MeshSnapshot(run.interior_material if len(run.interior_material[0]) else None)
                require(not boundary_material.explicit_empty and not interior_material.explicit_empty,
                        'NOT_SUPPORTED: dual material component geometry is empty')
                dual_material_audit = _dual_material_independent_audit(
                    boundary_material, interior_material, boundary_density, interior_density,
                    diagnostics['geometry'], source_settings.root_side)
            design = problem.design_regions[0]
            clearance = n.check_exact_boundary_clearance(design.mesh.exact_region(),cavity.exact_region(),0.)
            require(clearance['completed'] and (cavity.explicit_empty or clearance['contact'] is False),
                    'NOT_SUPPORTED: cavity touches original design boundary')
            retained = _append_reversed(design.mesh, cavity)
            require(diagnostics['solve']['status'] == 'BINARY_MOMENT_VALIDATED', 'geometry has no numerical success')
            claim = _normalize_model_moments(diagnostics['geometry']['source_subtraction_moments'],problem.normalization)
            artifact = SealedSpinArtifact.create(problem, outer=outer,
                design_results=[(design.region_id,retained,cavity)], final_material=material, newly_carved=cavity,
                solver_claim={'frame':FRAME,'order':list(ORDER),'moments':claim})
        except (ValueError, RuntimeError) as exc:
            error = str(exc)
    elif diagnostics['solve']['status'] == 'BINARY_MOMENT_VALIDATED':
        error = 'native geometry unavailable: '+diagnostics['geometry']['status']+': '+diagnostics['geometry']['diagnostic']
    identity = {'schema':'source-spin2014-locked-material-run-v1','problem_digest':problem.input_digest,
        'input_identity':problem.identity,'source_settings':asdict(source_settings),'solver_settings':asdict(solver_settings),
        'boundary_settings':asdict(boundary_settings),'density':boundary_density,
        'material_densities':{'boundary':boundary_density,'interior':interior_density},
        'minimum_outer_wall':problem.minimum_outer_wall,
        'settings_authority':'IMPLEMENTATION_DEFAULT','geometry_frame':'MODEL_COORDINATES_NATIVE_Z_ORIGIN_ZERO',
        'claim_frame':FRAME,'normalization':problem.normalization,'native_run':diagnostics,
        'status':diagnostics['solve']['status'],'seal_error':error,
        'artifact_scientific_digest':None if artifact is None else artifact.scientific_digest,
        'outer':None if outer is None else outer.identity,'material':None if material is None else material.identity,
        'cavity':None if cavity is None else cavity.identity,
        'boundary_material':None if boundary_material is None else boundary_material.identity,
        'interior_material':None if interior_material is None else interior_material.identity,
        'dual_material_independent_audit':dual_material_audit,
        'project_status':'PARTIAL_BLOCKED',
        'solver_target_feasibility':'UNKNOWN','physical_validation':'NOT_IN_SCOPE',
        'original_evidence_admission':'BLOCKED_UNKNOWN','build_source_attestation':'NOT_VERIFIED'}
    after = _observation()
    evidence = {'scope':'OBSERVED_FILES_ON_DISK_ONLY','build_source_attestation':'NOT_VERIFIED',
                'before':before,'after':after,'freshness':'CURRENT' if before==after else 'STALE_RESULT',
                'native_budget_scope':'NATIVE_SOURCE_MATERIAL_CHAIN_ONLY_EXCLUDES_INDEPENDENT_PYTHON_SEAL'}
    result = object.__new__(SourceSpinMaterialResult)
    for name,value in dict(status=identity['status'],artifact=artifact,seal_error=error,
        scientific_digest=digest(identity),_identity=canonical(identity),_outer=outer,_material=material,_cavity=cavity,
        _boundary_material=boundary_material,_interior_material=interior_material,
        _evidence=canonical(evidence),_checkpoint=getattr(run,'checkpoint',None)).items():
        object.__setattr__(result,name,value)
    return result
