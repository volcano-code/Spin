"""Read-only exact vertex-link evidence; never repairs or admits a mesh."""
import argparse
import hashlib
import json
from pathlib import Path


def link_components(faces, vertex):
    adjacency = {}
    incident = []
    for index, face in enumerate(faces):
        if vertex not in face:
            continue
        incident.append(index)
        others = [v for v in face if v != vertex]
        if len(others) != 2 or others[0] == others[1]:
            raise ValueError('DEGENERATE_VERTEX_STAR')
        a, b = others
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    unseen = set(adjacency)
    components = []
    while unseen:
        todo = [min(unseen)]
        found = set()
        while todo:
            current = todo.pop()
            if current in found:
                continue
            found.add(current)
            todo.extend(adjacency[current] - found)
        unseen -= found
        components.append(sorted(found))
    return dict(incident_face_indices=incident, components=components,
                link_degrees={str(v): len(adjacency[v]) for v in sorted(adjacency)},
                is_single_cycle=bool(adjacency) and len(components) == 1 and
                all(len(edges) == 2 for edges in adjacency.values()))


def main(packet_path, audit_path, output):
    packet = json.loads(packet_path.read_text())
    audit = json.loads(audit_path.read_text())
    witnesses = []
    for vertex in audit['non_manifold_vertices']:
        evidence = link_components(packet['faces'], vertex)
        indices = sorted({vertex} | {i for c in evidence['components'] for i in c})
        witnesses.append(dict(vertex=vertex, exact_coordinate=packet['vertices'][vertex],
            local_exact_coordinates={str(i): packet['vertices'][i] for i in indices},
            incident_faces=[packet['faces'][i] for i in evidence['incident_face_indices']],
            **evidence))
    result = dict(schema='exact-vertex-link-witness-v1',
        input_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
        audit_sha256=hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        witnesses=witnesses, geometry_admission=False, optimization_admission=False,
        repair_applied=False)
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps([dict(vertex=w['vertex'], components=len(w['components']),
        component_sizes=[len(c) for c in w['components']],
        single_cycle=w['is_single_cycle']) for w in witnesses]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('packet', 'audit', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.packet, args.audit, args.output)
