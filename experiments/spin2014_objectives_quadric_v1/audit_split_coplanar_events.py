"""Exact, saved-material-only witnesses for split queue 04 topology events.

No native mapping, Boolean or optimizer calls. Does not mutate the frozen
aggregator or grant its failed derivative gate. Run via source_job guardian.
"""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys
import gmpy2

from production_point_job import require_point_resource_envelope, save
from production_status import digest
from make_it_stand.verification.coplanar_patch import certify_convex_quad_retriangulation
from make_it_stand.verification.quadric_point_jobs import _read, verify_point_sources


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--job-envelope', type=Path, required=True)
    args = parser.parse_args()
    require_point_resource_envelope(args.job_envelope)
    sys.path.insert(0, str(args.evidence.resolve()))
    import diagnose_independent_primal_secants as primal
    import diagnose_independent_transport as transport
    import audit_production_coplanar_patches as patches
    manifest_path = args.evidence / 'split-production-manifest-04/manifest.json'
    manifest = _read(manifest_path)
    verify_point_sources(manifest, manifest['mapping_settings']['source_paths'])
    args.output.mkdir(exist_ok=False)
    indices = [0, *range(57, 65), *range(73, 81)]
    bound = {str(manifest_path.resolve()): digest(manifest_path),
        str(Path(__file__).resolve()): digest(__file__),
        str(Path(primal.__file__).resolve()): digest(primal.__file__),
        str(Path(transport.__file__).resolve()): digest(transport.__file__),
        str(Path(patches.__file__).resolve()): digest(patches.__file__)}

    def load(index):
        root = args.evidence / f'split-production-point-04-{index:03d}/output'
        receipt = _read(root / 'terminal.json')
        path = root / 'point.json'
        if receipt['status'] != 'COMPLETE' or receipt['index'] != index or \
                receipt['manifest_sha256'] != manifest['identity'] or digest(path) != receipt['packet_sha256']:
            raise ValueError('point packet is incomplete or changed')
        bound[str(path.resolve())] = digest(path)
        return _read(path)['evaluation']

    def coords(rows):
        return [[gmpy2.mpq(x) for x in row] for row in rows]

    base = load(0)
    df, sf = base['design_faces'], base['source_faces']
    design = coords(base['design_vertices_target_authority'])
    identities = primal.descriptors(design, df, coords(base['source_vertices_target_authority']),
        sf, coords(base['material_vertices_target_authority']),
        base['material_faces'], base['material_birth_faces'])
    if len(set(identities)) != len(identities):
        raise ValueError('base feature identities are not unique')
    metadata = ('surface_type', 'effective_normal_subdivisions', 'source_truncation_radius',
        'source_interface_face_count', 'source_cap_face_count', 'boolean_backend_identity')
    rows = []
    for index in indices:
        value = base if index == 0 else load(index)
        if value['source_faces'] != sf or value['design_vertices_target_authority'] != base['design_vertices_target_authority'] or \
                value['design_faces'] != df or any(value[k] != base[k] for k in metadata):
            raise ValueError('unexplained source/domain/construction metadata change')
        source = coords(value['source_vertices_target_authority'])
        points = coords(value['material_vertices_target_authority'])
        own = primal.descriptors(design, df, source, sf, points,
            value['material_faces'], value['material_birth_faces'])
        if len(set(own)) != len(own) or set(own) != set(identities):
            raise ValueError('minimum-dimensional features actually changed')
        if patches.topology_hash(own, value['material_faces'], value['material_birth_faces']) != value['provenance_topology_hash']:
            raise ValueError('saved provenance topology preimage mismatch')
        mapped = primal.reconstruct(design, df, source, sf, identities)
        lookup = {tuple(point): i for i, point in enumerate(points)}
        if len(lookup) != len(mapped) or set(map(tuple, mapped)) != set(lookup):
            raise ValueError('independently reconstructed intersections differ')
        order = [lookup[tuple(point)] for point in mapped]
        replay = [[order[i] for i in face] for face in base['material_faces']]
        expected = Counter((birth, patches.triangle(face)) for birth, face in zip(base['material_birth_faces'], replay))
        actual = Counter((birth, patches.triangle(face)) for birth, face in zip(value['material_birth_faces'], value['material_faces']))
        first, second = expected - actual, actual - expected
        proofs = []
        for birth in sorted({b for b, _ in first} | {b for b, _ in second}):
            if birth < len(df):
                raise ValueError('fixed-domain topology changed')
            f1 = [list(f) for (b, f), count in first.items() if b == birth for _ in range(count)]
            f2 = [list(f) for (b, f), count in second.items() if b == birth for _ in range(count)]
            corners = sorted({i for f in f1 + f2 for i in f})
            local = {old: new for new, old in enumerate(corners)}
            exact = [[Fraction(x.numerator, x.denominator) for x in points[i]] for i in corners]
            proof = certify_convex_quad_retriangulation(exact,
                [[local[i] for i in f] for f in f1], [[local[i] for i in f] for f in f2])
            for i in corners:
                transport.feature(source, sf, [birth - len(df)], points[i])
            proofs.append(dict(source_face=birth-len(df), corners=corners, proof=proof))
        row = dict(index=index, material_scientific_hash=value['material_scientific_hash'],
            provenance_topology_hash=value['provenance_topology_hash'],
            same_exact_vertex_features=True, independently_reconstructed_vertices=True,
            same_parameter_oriented_entity=True, patches=proofs)
        save(args.output / f'certificate-{index:02d}.json', row)
        rows.append(row)
        print(json.dumps(dict(index=index, certified=True, patches=len(proofs))), flush=True)
    verify_point_sources(manifest, manifest['mapping_settings']['source_paths'])
    if any(digest(path) != expected for path, expected in bound.items()):
        raise ValueError('proof inputs changed during inspection')
    save(args.output / 'terminal.json', dict(schema='split-production-coplanar-event-witness-v1',
        manifest_identity=manifest['identity'], sources=bound, rows=rows,
        native_invocations=0, boolean_invocations=0, optimizer_invocations=0,
        scientific_gate_granted=False, original_failed_aggregation_preserved=True))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as error:
        if '--output' in sys.argv:
            output = Path(sys.argv[sys.argv.index('--output') + 1])
            if output.is_dir() and not (output / 'failure.json').exists():
                save(output / 'failure.json', dict(schema='split-production-coplanar-event-failure-v1',
                    exception_type=type(error).__name__, exception=str(error),
                    scientific_gate_granted=False, incomplete_prefix_preserved=True))
        raise
