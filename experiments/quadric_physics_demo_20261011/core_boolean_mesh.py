"""Zero-displacement binary64 vertex canonicalization for core Boolean output."""
import numpy as np
import re


def binary64_rational_token(value):
    """Exact finite binary64 value, including native-required /1 for integers."""
    value=float(value)
    if not np.isfinite(value):raise ValueError('NONFINITE_RATIONAL_COORDINATE')
    numerator,denominator=value.as_integer_ratio()
    return str(numerator)+'/'+str(denominator)


def manifest_design_seed(manifest):
    """Select the declared domain operand; the outer surface is not a fallback."""
    expression=manifest.get('construction',{}).get('design_operation','')
    match=re.fullmatch(r'DIFFERENCE\(([A-Za-z0-9_-]+), UNION\(ALL fixed_parts\)\)',expression)
    if match is None:raise ValueError('UNSUPPORTED_MANIFEST_DESIGN_OPERATION')
    name=match.group(1)
    if name not in manifest.get('parts',{}):raise ValueError('MISSING_MANIFEST_DESIGN_SEED')
    return name
import trimesh


def bounded_rational_core_mesh(rational_vertices, faces, *, length_scale,
                              boundary_witness_policy='SURVIVING_VERTEX_V1'):
    """Cancel only opposite faces after identical binary64-coordinate collapse.

    Exact L1 bounds certify the coordinate conversion and every erased patch
    against a surviving vertex. The swept-volume bound covers triangle tubes
    and erased-patch boxes. Set validity/self-intersection still require the
    independent checker; this function never grants material admission.
    """
    import gmpy2
    if boundary_witness_policy not in ('SURVIVING_VERTEX_V1','SURVIVING_EDGE_V2'):
        raise ValueError('UNKNOWN_BOUNDARY_WITNESS_POLICY')
    if not np.isfinite(length_scale) or length_scale <= 0:
        raise ValueError('INVALID_FROZEN_LENGTH_SCALE')
    exact = [tuple(gmpy2.mpq(x) for x in row) for row in rational_vertices]
    vertices = np.asarray([[float(x) for x in row] for row in exact])
    faces = np.asarray(faces)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or not np.isfinite(vertices).all():
        raise ValueError('INVALID_CORE_VERTICES')
    if faces.ndim != 2 or faces.shape[1] != 3 or faces.dtype.kind not in 'iu' or not len(faces):
        raise ValueError('INVALID_CORE_FACES')
    if faces.min() < 0 or faces.max() >= len(vertices):
        raise ValueError('CORE_FACE_INDEX_OUT_OF_RANGE')
    represented = [tuple(gmpy2.mpq(float(x)) for x in row) for row in vertices]
    delta = max(sum(abs(a-b) for a,b in zip(p,q)) for p,q in zip(exact,represented))
    unique, inverse = np.unique(vertices, axis=0, return_inverse=True)
    mapped = inverse[faces]
    collapsed = np.any(np.diff(np.sort(mapped, axis=1), axis=1) == 0, axis=1)
    groups = {}
    for index in np.flatnonzero(~collapsed):
        triangle = mapped[index]
        key = tuple(sorted(int(x) for x in triangle))
        parity = sum(int(triangle[a] > triangle[b]) for a,b in ((0,1),(0,2),(1,2))) % 2
        groups.setdefault(key, []).append((int(index), parity))
    removed = set(int(x) for x in np.flatnonzero(collapsed))
    pairs = []
    for occurrences in groups.values():
        if len(occurrences) == 1:
            continue
        if len(occurrences) != 2 or occurrences[0][1] == occurrences[1][1]:
            raise ValueError('NONCANCELLING_DUPLICATE_CORE_FACES')
        pair = [x[0] for x in occurrences]
        pairs.append(pair); removed.update(pair)
    keep = [i for i in range(len(faces)) if i not in removed]
    if not keep:
        raise ValueError('EMPTY_CORE_BOUNDARY')
    used = set(int(x) for x in mapped[keep].ravel())
    surviving_edges=set()
    if boundary_witness_policy=='SURVIVING_EDGE_V2':
        for triangle in mapped[keep]:
            surviving_edges.update(tuple(sorted((int(triangle[i]),int(triangle[(i+1)%3])))) for i in range(3))
    witnesses = []
    boundary_error = delta
    patch_volume = gmpy2.mpq(0)
    for index in sorted(removed):
        candidates = sorted(set(int(x) for x in mapped[index]) & used)
        if not candidates:
            raise ValueError('ERASED_PATCH_HAS_NO_SURVIVING_BOUNDARY_WITNESS')
        best = None
        for candidate in candidates:
            point = tuple(gmpy2.mpq(float(x)) for x in unique[candidate])
            radius = max(sum(abs(a-b) for a,b in zip(exact[int(vertex)],point))
                         for vertex in faces[index])
            if best is None or radius < best[0]:
                best = (radius,candidate)
        radius,candidate = best
        witness=dict(face=index,surviving_coordinate=unique[candidate].tolist(),
                     distance_l1_upper_bound_exact=str(radius))
        erased_volume=8*radius**3
        if boundary_witness_policy=='SURVIVING_EDGE_V2':
            points=[exact[int(vertex)] for vertex in faces[index]]
            for a in candidates:
                for b in candidates:
                    if a>=b or (a,b) not in surviving_edges:continue
                    start=tuple(gmpy2.mpq(float(x)) for x in unique[a])
                    end=tuple(gmpy2.mpq(float(x)) for x in unique[b])
                    direction=tuple(y-x for x,y in zip(start,end))
                    denominator=sum(x*x for x in direction)
                    edge_radius=gmpy2.mpq(0)
                    for point in points:
                        t=min(gmpy2.mpq(1),max(gmpy2.mpq(0),
                            sum((x-y)*d for x,y,d in zip(point,start,direction))/denominator))
                        edge_radius=max(edge_radius,sum(abs(x-(s+t*d)) for x,s,d in zip(point,start,direction)))
                    if edge_radius>=radius:continue
                    radius=edge_radius
                    # Convexity bounds the whole source triangle by this capsule.
                    # Its reverse edge lies within coordinate-rounding delta of
                    # the original triangle. 8*r²*(L1 length+2r) bounds its volume.
                    edge_span=sum(abs(x) for x in direction)
                    erased_volume=8*radius**2*(edge_span+2*radius)
                    witness=dict(face=index,surviving_edge=[unique[a].tolist(),unique[b].tolist()],
                        surviving_edge_indices=[a,b],distance_l1_upper_bound_exact=str(radius),
                        capsule_volume_upper_bound_exact=str(erased_volume))
        boundary_error = max(boundary_error,radius)
        patch_volume += erased_volume
        witnesses.append(witness)
    # L1 span bounds every original triangle edge and its area/perimeter.
    span = sum(max(p[a] for p in exact)-min(p[a] for p in exact) for a in range(3))
    swept = len(faces)*(delta*span**2 + 12*delta**2*span + 6*delta**3)+patch_volume
    length = gmpy2.mpq(float(length_scale))
    if boundary_error > length/gmpy2.mpq(10**6) or swept > length**3/gmpy2.mpq(10**6):
        raise ValueError('BOUNDED_CORE_REPRESENTATION_ERROR_EXCEEDED')
    mesh, report = canonical_core_mesh(vertices,faces[keep])
    def upper(value):
        number = float(value)
        return float(np.nextafter(number,np.inf)) if gmpy2.mpq(number) < value else number
    report.update(profile='BOUNDED_RATIONAL_ALIAS_CHAIN_V1',
        removed_collapsed_faces=int(np.count_nonzero(collapsed)),
        removed_opposite_face_pairs=len(pairs), opposite_face_pairs=pairs,
        retained_original_face_indices=keep, removed_face_witnesses=witnesses,
        coordinate_rounding_l1_upper_bound_exact=str(delta),
        maximum_vertex_displacement=upper(delta),
        maximum_vertex_displacement_frame='rational_authority_to_binary64',
        boundary_deviation_upper_bound_exact=str(boundary_error),
        boundary_deviation_upper_bound=upper(boundary_error),
        swept_volume_upper_bound_exact=str(swept), swept_volume_upper_bound=upper(swept),
        frozen_length_scale=length_scale, boundary_correspondence_certified=True,
        volume_bound_scope='swept_closed_chain_bound; region validity requires independent audit',
        independent_material_admission=False)
    if boundary_witness_policy=='SURVIVING_EDGE_V2':
        report.update(profile='BOUNDED_RATIONAL_ALIAS_EDGE_WITNESS_V2',
                      boundary_witness_policy=boundary_witness_policy)
    return mesh,report


def canonical_core_mesh(vertices,faces):
    vertices=np.asarray(vertices,dtype=np.float64)
    faces=np.asarray(faces)
    if vertices.ndim!=2 or vertices.shape[1]!=3 or not np.isfinite(vertices).all():
        raise ValueError('INVALID_CORE_VERTICES')
    if faces.ndim!=2 or faces.shape[1]!=3 or faces.dtype.kind not in 'iu' or not len(faces):
        raise ValueError('INVALID_CORE_FACES')
    if faces.min()<0 or faces.max()>=len(vertices):raise ValueError('CORE_FACE_INDEX_OUT_OF_RANGE')
    unique,inverse=np.unique(vertices,axis=0,return_inverse=True)
    remapped=inverse[faces]
    collapsed=(remapped[:,0]==remapped[:,1]) | (remapped[:,1]==remapped[:,2]) | (remapped[:,2]==remapped[:,0])
    retained=remapped[~collapsed]
    if not len(retained):raise ValueError('EMPTY_CORE_BOUNDARY')
    used=np.unique(retained)
    compact=np.full(len(unique),-1,dtype=np.int64);compact[used]=np.arange(len(used))
    mesh=trimesh.Trimesh(unique[used],compact[retained],process=False)
    if np.any(mesh.area_faces==0):raise ValueError('COLLINEAR_DISTINCT_VERTEX_FACE')
    if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume<=0:
        raise ValueError('CORE_BOUNDARY_NOT_CLOSED_ORIENTED_POSITIVE')
    # Collapsing coordinate aliases must not turn two vertex fans into a pinch.
    counts=np.bincount(inverse,minlength=len(unique))
    affected=compact[np.flatnonzero(counts>1)];affected=affected[affected>=0]
    for vertex in affected:
        incident=mesh.faces[np.any(mesh.faces==vertex,axis=1)]
        links={}
        for triangle in incident:
            a,b=[int(x) for x in triangle if x!=vertex]
            links.setdefault(a,[]).append(b);links.setdefault(b,[]).append(a)
        if not links or any(len(neighbors)!=2 for neighbors in links.values()):
            raise ValueError('CORE_COLLAPSE_NONMANIFOLD_VERTEX_LINK')
        reached=set();todo=[next(iter(links))]
        while todo:
            current=todo.pop()
            if current in reached:continue
            reached.add(current);todo.extend(links[current])
        if len(reached)!=len(links):raise ValueError('CORE_COLLAPSE_DISCONNECTED_VERTEX_LINK')
    return mesh,dict(profile='IDENTICAL_BINARY64_VERTEX_CANONICALIZATION_V1',
        merged_coordinate_aliases=int(len(vertices)-len(unique)),
        removed_collapsed_faces=int(np.count_nonzero(collapsed)),
        maximum_vertex_displacement=0.,independent_material_admission=False)
