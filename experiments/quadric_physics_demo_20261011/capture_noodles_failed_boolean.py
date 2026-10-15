"""Capture the frozen 73rd Boolean for diagnosis; never admit or repair it."""
import argparse
import hashlib
import json
from pathlib import Path


def main(case, predecessor, output):
    import gmpy2
    import numpy as np
    import trimesh
    import _make_it_stand as native
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    def save(name, data):
        def array_value(value):
            if isinstance(value, np.ndarray): return value.tolist()
            if isinstance(value, np.generic): return value.item()
            raise TypeError('UNSUPPORTED_DIAGNOSTIC_VALUE:'+type(value).__name__)
        with (output/name).open('x') as stream:
            json.dump(data, stream, default=array_value, allow_nan=False)
    manifest = json.loads((case/'manifest.json').read_text())
    packet = json.loads(predecessor.read_text())
    prefix = packet['excluded_parts']
    if len(prefix) != 72 or prefix != manifest['construction']['fixed_parts'][:72] or packet['manifest_sha256'] != sha(case/'manifest.json'):
        raise ValueError('FROZEN_FAILED_BOOLEAN_PREDECESSOR_MISMATCH')
    name = manifest['construction']['fixed_parts'][72]
    entry = manifest['parts'][name]
    path = case/entry['model_obj']
    if sha(path) != entry['sha256']['model']:
        raise ValueError('FIXED_PART_CHANGED')
    output.mkdir(parents=True, exist_ok=False)
    part = trimesh.load(path, force='mesh', process=False)
    result = native.libigl_exact_boolean_construction(np.asarray(packet['vertices']),
        np.asarray(packet['faces'], dtype=np.int32), part.vertices, part.faces, 'difference')
    # Persist before any transformation or validity claim.
    save('original-rational-output.json', dict(result))
    if not result['succeeded']:
        raise ValueError('ACTUAL_BOOLEAN_FAILED')
    vertices = np.asarray([[float(gmpy2.mpq(x)) for x in row] for row in result['vertices']])
    faces = np.asarray(result['faces'], dtype=np.int32)
    save('binary64-output.json', dict(vertices=vertices.tolist(), faces=faces.tolist()))
    mesh = trimesh.Trimesh(vertices, faces, process=False)
    unique, inverse = np.unique(vertices, axis=0, return_inverse=True)
    mapped = inverse[faces]
    collapsed = np.any(np.diff(np.sort(mapped, axis=1), axis=1) == 0, axis=1)
    retained = mapped[~collapsed]
    groups = {}
    for index, triangle in enumerate(retained):
        key = tuple(sorted(int(x) for x in triangle))
        parity = sum(int(triangle[a] > triangle[b]) for a, b in ((0, 1), (0, 2), (1, 2))) % 2
        groups.setdefault(key, []).append(dict(retained_index=index, parity=parity))
    edges = np.sort(np.concatenate((retained[:,[0,1]], retained[:,[1,2]], retained[:,[2,0]])), axis=1)
    edge_keys, counts = np.unique(edges, axis=0, return_counts=True)
    save('diagnosis.json', dict(part=name, predecessor_sha256=sha(predecessor),
        native_sha256=sha(Path(native.__file__)), zero_area_face_indices=np.flatnonzero(mesh.area_faces == 0).tolist(),
        coordinate_alias_count=int(len(vertices)-len(unique)), collapsed_face_indices=np.flatnonzero(collapsed).tolist(),
        duplicate_triangle_groups=[dict(vertices=list(key), occurrences=value) for key, value in groups.items() if len(value)>1],
        abnormal_edges=[dict(vertices=key.tolist(), incidence=int(count)) for key, count in zip(edge_keys,counts) if count != 2],
        geometry_admission=False, optimization_admission=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--predecessor', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.case, args.predecessor, args.output)
