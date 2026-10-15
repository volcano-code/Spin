"""Derived binary64 exporter: use the existing complete material gate unchanged.

Exact authority is preserved. Only a new binary64 derivative is proposed;
moments alone or local patch proofs never qualify it as a material object.
The caller persists both objects independently and binds this attachment.
"""
from fractions import Fraction as Q
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.spin_artifact import SealedSpinArtifact
from make_it_stand.verification.spin_problem import MeshSnapshot, partition
from make_it_stand.verification.spin_common import digest
from make_it_stand.verification.rational_text import fraction_from_rational_text
from ._coplanar_boundary import propose_coplanar_diagonal_repairs


def derive_float_material(authority):
    if not isinstance(authority,SealedSpinArtifact):
        raise TypeError('expected an actual sealed authority')
    authority=authority.revalidate()
    native=spin._native_module()
    records=[];cache={}

    def derive(mesh,role):
        key=digest(mesh.identity)
        if key in cache:
            candidate,source_role=cache[key]
            records.append({'role':role,'reuses_identical_source_role':source_role,
                            'authority_mesh_digest':key})
            return candidate
        if mesh.explicit_empty:
            candidate=MeshSnapshot(None)
            records.append({'role':role,'explicit_empty':True,'authority_mesh_digest':key})
            cache[key]=(candidate,role)
            return candidate
        vertices,faces=mesh.arrays()
        before=dict(native.audit_exact_triangle_boundary(vertices,faces,'AnyNonzero'))
        exact_rows=mesh.rational_vertices
        exact=([[fraction_from_rational_text(x) for x in row] for row in exact_rows]
               if exact_rows is not None else [[Q(float(x)) for x in row] for row in vertices])
        hints=list(before['intersecting_pairs'])
        # The native checker correctly skips intersection predicates if any
        # represented face is degenerate. Its actual degenerate-face list is
        # therefore a separate source of hints, never proof of a safe repair.
        if before['degenerate_faces']:
            adjacency={}
            for index,face in enumerate(faces.tolist()):
                for a,b in zip(face,face[1:]+face[:1]):
                    adjacency.setdefault(tuple(sorted((a,b))),[]).append(index)
            for index in before['degenerate_faces']:
                face=faces[index].tolist()
                for a,b in zip(face,face[1:]+face[:1]):
                    hints.extend([index,other] for other in adjacency[tuple(sorted((a,b)))] if other!=index)
        proposal=propose_coplanar_diagonal_repairs(exact,vertices.tolist(),faces.tolist(),
            hints,face_regions=[role]*len(faces))
        # MeshSnapshot rechecks nesting, winding and signed components via Nef;
        # this is stronger than accepting the AnyNonzero boundary report above.
        candidate=MeshSnapshot((vertices,np.asarray(proposal['faces'],dtype=np.int32)))
        error2=max(sum((Q(float(x))-e)**2 for x,e in zip(v,row))
                   for v,row in zip(vertices,exact))
        records.append({'role':role,'authority_mesh_digest':key,
            'derived_mesh_digest':digest(candidate.identity),'before_boundary':before,
            'proposal':proposal,'maximum_corresponding_vertex_error_squared_exact':str(error2),
            'relation':'EXACT_PATCH_RETRIANGULATION_THEN_RECORDED_BINARY64_ROUNDING',
            'source_face_provenance':'NOT_USED_FOR_OUTPUT_ONLY_EXACT_PATCH_PROOF'})
        cache[key]=(candidate,role)
        return candidate

    # Fixed and outer input identities are immutable; do not silently replace
    # these parts of the optimization problem with rounded/modified operands.
    for region in authority.problem.fixed_regions:
        vertices,faces=region.mesh.arrays()
        if region.mesh.explicit_empty:continue
        facts=dict(native.audit_exact_triangle_boundary(vertices,faces,'AnyNonzero'))
        if not facts['boundary_valid']:
            raise ValueError('fixed input has no valid binary64 boundary: '+region.region_id)
        if region.mesh.rational_vertices is not None:
            # A valid display alone is not proof that the unchanged fixed
            # material represents the authority used in the partition check.
            partition(region.mesh.exact_region(),[MeshSnapshot((vertices,faces)).exact_region()])
    outer_vertices,outer_faces=authority.outer.arrays()
    if not native.audit_exact_triangle_boundary(outer_vertices,outer_faces,'AnyNonzero')['boundary_valid']:
        raise ValueError('immutable outer input has no valid binary64 boundary')
    if authority.outer.rational_vertices is not None:
        partition(authority.outer.exact_region(),
                  [MeshSnapshot((outer_vertices,outer_faces)).exact_region()])
    results=[(name,derive(retained,'design.'+name+'.retained'),
              derive(carved,'design.'+name+'.carved'))
             for name,retained,carved in authority.design_results]
    material=derive(authority.final_material,'final.material')
    carved=derive(authority.newly_carved,'final.newly_carved')
    # This existing gate checks exact partitions, containment, signed regions,
    # clearance and independent density-aware integrals. It can reject a tiny
    # rounding gap even when all separate boundaries and mass checks pass.
    candidate=SealedSpinArtifact.create(authority.problem,outer=authority.outer,
        design_results=results,final_material=material,newly_carved=carved,
        solver_claim=authority.solver_claim)
    attachment={'schema':'derived-coplanar-material-v1',
        'scope':'DERIVED_MESH_MATERIAL_NOT_MANUFACTURING_OR_SOLVER_APPROVAL',
        'authority_scientific_digest':authority.scientific_digest,
        'derived_scientific_digest':candidate.scientific_digest,
        'problem_digest':authority.problem.input_digest,'records':records,
        'manufacturing_admission':candidate.facts['manufacturing_admission']}
    return candidate,attachment
