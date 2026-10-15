"""SYNTHETIC_REGRESSION for the pinned Cherchi Boolean component only."""

from __future__ import annotations

import hashlib
from fractions import Fraction
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from make_it_stand.verification.mass_moments import TriangleMesh, integrate_geometric_moments
from make_it_stand.verification.cherchi_face_provenance import recover_face_sources


def exact_binary_rational(value: str) -> Fraction:
    if value == "0/1":
        return Fraction(0)
    negative = value.startswith("-")
    if negative:
        value = value[1:]
    numerator, denominator = value.split("/")

    def integer(token: str) -> int:
        return int(token[2:], 2) if token.startswith("0b") else int(token)

    result = Fraction(integer(numerator), integer(denominator))
    return -result if negative else result


def box(path: Path, lo: tuple[float, ...], hi: tuple[float, ...]) -> None:
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    vertices = ((x0,y0,z0), (x1,y0,z0), (x1,y1,z0), (x0,y1,z0),
                (x0,y0,z1), (x1,y0,z1), (x1,y1,z1), (x0,y1,z1))
    faces = ((0,2,1),(0,3,2), (4,5,6),(4,6,7),
             (0,1,5),(0,5,4), (1,2,6),(1,6,5),
             (2,3,7),(2,7,6), (3,0,4),(3,4,7))
    with path.open("w", encoding="ascii") as stream:
        for vertex in vertices:
            stream.write("v " + " ".join(map(str, vertex)) + "\n")
        for face in faces:
            stream.write("f " + " ".join(str(index + 1) for index in face) + "\n")


def read_obj(path: Path) -> TriangleMesh:
    vertices = []
    faces = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "v":
            vertices.append(tuple(float(value) for value in fields[1:4]))
        elif fields[0] == "f":
            faces.append(tuple(int(value.split("/")[0]) - 1 for value in fields[1:4]))
    return TriangleMesh(tuple(vertices), tuple(faces))


def write_mesh(path: Path, vertices, faces) -> None:
    with path.open("w", encoding="ascii") as stream:
        for vertex in vertices:
            stream.write("v " + " ".join(format(value, ".17g") for value in vertex) + "\n")
        for face in faces:
            stream.write("f " + " ".join(str(index + 1) for index in face) + "\n")


def main() -> None:
    executable, work = Path(sys.argv[1]), Path(sys.argv[2])
    work.mkdir(parents=True, exist_ok=False)
    left, right = work / "left.obj", work / "right.obj"
    box(left, (0, 0, 0), (3, 2, 1.5))
    box(right, (1, .5, .25), (4, 2.5, 1.25))
    assert abs(integrate_geometric_moments(read_obj(left)).volume - 9) < 1e-12
    assert abs(integrate_geometric_moments(read_obj(right)).volume - 6) < 1e-12
    digests = []
    measured = []
    for repetition in (1, 2):
        output = work / f"intersection-{repetition}.obj"
        result = subprocess.run((str(executable), "intersection", str(left), str(right), str(output)),
                                check=True, capture_output=True, text=True, timeout=120)
        print(f"run={repetition} stdout={result.stdout.strip()!r}")
        digests.append(hashlib.sha256(output.read_bytes()).hexdigest())
        moments = integrate_geometric_moments(read_obj(output))
        expected = (3.0, (6.0, 3.75, 2.25),
                    ((13.0, 7.5, 4.5), (7.5, 5.25, 2.8125), (4.5, 2.8125, 1.9375)))
        observed = (moments.volume, moments.first_moment, moments.second_moment)
        measured.append(observed)
        for actual, target in zip((observed[0], *observed[1], *(v for row in observed[2] for v in row)),
                                  (expected[0], *expected[1], *(v for row in expected[2] for v in row))):
            assert abs(actual - target) <= 1e-9 * max(1.0, abs(target)), (actual, target)
        print(f"run={repetition} vertices={len(read_obj(output).vertices)} faces={len(read_obj(output).faces)} "
              f"volume={moments.volume:.17g} first_moment={moments.first_moment} sha256={digests[-1]}")
    assert measured[0] == measured[1], measured
    assert digests[0] == digests[1], digests

    # The exact diagnostic path exports nearest binary64 coordinates of the
    # rational construction points. Its bytes may intentionally differ from
    # the legacy approximate export, but must be deterministic and auditable.
    exact_output = work / "intersection-exact.obj"
    exact_run = subprocess.run(
        (str(executable), "intersection", str(left), str(right), str(exact_output),
         "--dump-exact-rational"),
        check=True, capture_output=True, text=True, timeout=120,
    )
    exact_repeat_output = work / "intersection-exact-repeat.obj"
    subprocess.run(
        (str(executable), "intersection", str(left), str(right),
         str(exact_repeat_output), "--dump-exact-rational"),
        check=True, capture_output=True, text=True, timeout=120,
    )
    assert exact_output.read_bytes() == exact_repeat_output.read_bytes()
    assert exact_repeat_output.read_bytes() == exact_output.read_bytes()
    for suffix in (".exact.tsv", ".construction.tsv", ".faces.tsv"):
        assert Path(str(exact_repeat_output) + suffix).read_bytes() == \
               Path(str(exact_output) + suffix).read_bytes(), suffix
    lines = Path(str(exact_output) + ".exact.tsv").read_text(encoding="ascii").splitlines()
    assert lines[0].split("\t") == [
        "vertex_index", "status", "exact_x", "exact_y", "exact_z",
        "memory_x", "memory_y", "memory_z", "obj_x", "obj_y", "obj_z",
    ]
    exact_rows = [line.split("\t") for line in lines[1:]]
    exact_obj_vertices = read_obj(exact_output).vertices
    assert len(exact_rows) == len(exact_obj_vertices)
    for index, row in enumerate(exact_rows):
        assert int(row[0]) == index and row[1] in {
            "MATCHED", "APPROXIMATION_DISAGREES",
        }, row
        memory = tuple(float(value) for value in row[5:8])
        assert memory == exact_obj_vertices[index], (index, memory, exact_obj_vertices[index])
        if row[1] == "MATCHED":
            exact = tuple(float(exact_binary_rational(value)) for value in row[2:5])
            assert exact == memory, (index, exact, memory)
        else:
            assert any(
                float(exact_binary_rational(value)) != memory[axis]
                for axis, value in enumerate(row[2:5])
            ), row

    # The approximate construction-coordinate path can move a true Boolean
    # intersection just outside a design face after rescaling. The exact
    # diagnostic export must round the exact rational intersection to the
    # nearest binary64 coordinate, keeping the represented result in-domain.
    scale = 100.3
    domain_vertices = tuple(
        (scale*x, scale*y, scale*z)
        for x,y,z in ((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                       (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
    )
    domain_faces = ((0,2,1),(0,3,2),(4,5,6),(4,6,7),
                    (0,1,5),(0,5,4),(1,2,6),(1,6,5),
                    (2,3,7),(2,7,6),(3,0,4),(3,4,7))
    cutter_vertices = (
        (-38.689039912985244, 1.8093672852342975, -15.256169621957516),
        (196.4271841245415, -31.564373718853727, -3.755720797629939),
        (218.90931018006043, 0.7521466215747927, -35.01066273073553),
        (248.50539822054958, -40.038742338155544, 65.75607915778316),
    )
    cutter_faces = ((0,2,1),(0,1,3),(0,3,2),(1,2,3))
    domain_path, cutter_path = work/"scaled-domain.obj", work/"oblique-cutter.obj"
    rounded_output = work/"scaled-subtraction-exact.obj"
    write_mesh(domain_path, domain_vertices, domain_faces)
    write_mesh(cutter_path, cutter_vertices, cutter_faces)
    subprocess.run(
        (str(executable), "subtraction", str(domain_path), str(cutter_path),
         str(rounded_output), "--dump-exact-rational"),
        check=True, capture_output=True, text=True, timeout=120,
    )
    rounded_vertices = read_obj(rounded_output).vertices
    # Provenance for the *actual* Boolean mesh must not silently substitute
    # nearest-exact coordinates. The scaled oblique cut contains implicit
    # construction points for which those two exports can differ.
    ordinary_output = work / "scaled-subtraction-ordinary.obj"
    actual_provenance_output = work / "scaled-subtraction-provenance.obj"
    subprocess.run(
        (str(executable), "subtraction", str(domain_path), str(cutter_path),
         str(ordinary_output)),
        check=True, capture_output=True, text=True, timeout=120,
    )
    subprocess.run(
        (str(executable), "subtraction", str(domain_path), str(cutter_path),
         str(actual_provenance_output), "--dump-actual-provenance"),
        check=True, capture_output=True, text=True, timeout=120,
    )
    assert actual_provenance_output.read_bytes() == ordinary_output.read_bytes(), (
        "provenance diagnostics changed the represented Boolean result"
    )
    assert Path(str(actual_provenance_output) + ".construction.tsv").is_file()
    assert Path(str(actual_provenance_output) + ".faces.tsv").is_file()
    for vertex in rounded_vertices:
        assert all(-scale <= coordinate <= scale for coordinate in vertex), (
            vertex, scale,
        )
    boundary_rows = [line.split("\t") for line in
                     Path(str(rounded_output) + ".exact.tsv").read_text(encoding="ascii").splitlines()[1:]]
    assert len(boundary_rows) == len(rounded_vertices)
    exact_scale = Fraction.from_float(scale)
    checked_boundary_vertices = 0
    for index, row in enumerate(boundary_rows):
        assert int(row[0]) == index
        vertex = rounded_vertices[index]
        for axis, coordinate in enumerate(vertex):
            # The explicit exact-export mode must use nearest coordinates even
            # when an ill-conditioned approximate construction is farther away.
            # This is a conversion contract, NOT a topology/manufacturing pass.
            expected_coordinate = float(exact_binary_rational(row[2 + axis]))
            assert coordinate == expected_coordinate, (
                'exact export retained an approximate construction coordinate',
                index, axis, coordinate, expected_coordinate,
            )
            if abs(coordinate) == scale:
                assert row[1] == "MATCHED", (index, row[1])
                assert exact_binary_rational(row[2 + axis]) == \
                    (exact_scale if coordinate > 0 else -exact_scale), (index, row, coordinate)
                checked_boundary_vertices += 1
    assert checked_boundary_vertices > 0

    # The exact diagnostic must retain the construction relation for each
    # exported Cherchi vertex, not merely its rounded coordinate.  This is
    # evidence for a future fixed-topology differential; it is not itself a
    # derivative and unsupported constructions must remain explicit.
    construction_path = Path(str(exact_output) + ".construction.tsv")
    construction_lines = construction_path.read_text(encoding="ascii").splitlines()
    construction_header = construction_lines[0].split("\t")
    assert construction_header[:4] == [
        "vertex_index", "point_type", "operand_count", "parameter_t",
    ]
    construction_rows = [line.split("\t") for line in construction_lines[1:]]
    assert len(construction_rows) == len(exact_rows)
    assert len(construction_header) == 31
    expected_operand_counts = {
        "EXPLICIT3D": 1,
        "LPI": 5,
        "TPI": 9,
        "LNC": 2,
    }
    source_vertices = set()
    source_facets = set()
    for input_mesh in (left, right):
        input_vertices = []
        for line in input_mesh.read_text(encoding="ascii").splitlines():
            fields = line.split()
            if fields and fields[0] == "v":
                point = tuple(Fraction.from_float(float(value)) for value in fields[1:4])
                input_vertices.append(point)
                source_vertices.add(point)
            elif fields and fields[0] == "f":
                face = tuple(int(value.split("/")[0]) - 1 for value in fields[1:4])
                source_facets.add(tuple(sorted(input_vertices[vertex] for vertex in face)))
    for index, row in enumerate(construction_rows):
        assert int(row[0]) == index, row
        assert row[1] in expected_operand_counts, row
        operand_count = expected_operand_counts[row[1]]
        assert int(row[2]) == operand_count, row
        assert len(row) == len(construction_header), row
        if row[1] == "LNC":
            assert 0.0 <= float(row[3]) <= 1.0, row
        for coordinate_index, coordinate in enumerate(row[4:]):
            if coordinate_index >= 3 * operand_count:
                assert coordinate == "-", row
                continue
            exact_binary_rational(coordinate)
        if row[1] in {"LPI", "TPI", "LNC"}:
            for operand in range(operand_count):
                point = tuple(exact_binary_rational(row[4 + 3 * operand + axis])
                              for axis in range(3))
                assert point in source_vertices, (index, row[1], operand, point)
        if row[1] == "LPI":
            plane_points = tuple(sorted(
                tuple(exact_binary_rational(row[4 + 3 * operand + axis])
                      for axis in range(3))
                for operand in (2, 3, 4)
            ))
            assert plane_points in source_facets, (index, plane_points)
    face_lines = Path(str(exact_output) + ".faces.tsv").read_text(encoding="ascii").splitlines()
    assert face_lines[0].split("\t") == ["face_index", "input_object_mask"]
    face_rows = [line.split("\t") for line in face_lines[1:]]
    assert len(face_rows) == len(read_obj(exact_output).faces)
    object_masks = set()
    for index, row in enumerate(face_rows):
        assert int(row[0]) == index and row[1].startswith("0x"), row
        mask = int(row[1], 16)
        assert mask in {0x1, 0x2}, (index, hex(mask))
        object_masks.add(mask)
    assert object_masks == {0x1, 0x2}, object_masks
    recovered = recover_face_sources((left, right), exact_output, exact_rows, face_rows)
    assert len(recovered) == len(face_rows)
    for index, record in enumerate(recovered):
        assert record["face_index"] == index and record["source_triangles"], record
        label_mask = int(face_rows[index][1], 16)
        candidate_mask = 0
        for source in record["source_triangles"]:
            candidate_mask |= 1 << source["input_object"]
        assert candidate_mask & label_mask, (index, label_mask, record)
    print(f"exact_coordinate_diagnostic={exact_run.stdout.strip()} "
          f"sidecar_vertices={len(exact_rows)} sidecar_faces={len(face_rows)} "
          f"source_object_masks={sorted(hex(mask) for mask in object_masks)} "
          f"face_source_records={len(recovered)} "
          f"ambiguous_face_sources={sum(len(record['source_triangles']) != 1 for record in recovered)} "
          f"exact_export_repeatable=True")

    roundtrip_dir = work / "roundtrip-precision"
    roundtrip_dir.mkdir()
    precise_left, precise_right = roundtrip_dir / "left.obj", roundtrip_dir / "right.obj"
    box(precise_left, (0.12345678901234566, 0, 0), (3, 2, 1.5))
    box(precise_right, (1.2345678901234567, .5, .25), (4, 2.5, 1.25))
    precise_output = roundtrip_dir / "intersection.obj"
    subprocess.run(
        (str(executable), "intersection", str(precise_left), str(precise_right),
         str(precise_output), "--dump-exact-rational"),
        check=True, capture_output=True, text=True, timeout=120,
    )
    precise_rows = [line.split("\t") for line in
                    Path(str(precise_output) + ".exact.tsv").read_text(encoding="ascii").splitlines()[1:]]
    for row in precise_rows:
        assert row[1] == "MATCHED", row
        assert tuple(map(float, row[5:8])) == tuple(map(float, row[8:11])), row

    captured = Path(__file__).resolve().parent / 'fixtures/large-cap-cut-v1'
    captured_output = work / 'captured-large-cap-exact.obj'
    subprocess.run((str(executable), 'subtraction', str(captured/'design.obj'),
                    str(captured/'source.obj'), str(captured_output), '--dump-exact-rational'),
                   check=True, capture_output=True, text=True, timeout=120)
    captured_vertices = read_obj(captured_output).vertices
    captured_rows = [line.split('\t') for line in
                     Path(str(captured_output)+'.exact.tsv').read_text().splitlines()[1:]]
    assert len(captured_rows) == len(captured_vertices) > 0
    for index, row in enumerate(captured_rows):
        assert int(row[0]) == index
        expected = tuple(float(exact_binary_rational(token)) for token in row[2:5])
        assert tuple(captured_vertices[index]) == expected, (
            'captured exact export retained approximate coordinates', index,
            tuple(captured_vertices[index]), expected,
        )

    # A minus B: analytic raw moments of A minus analytic raw moments of overlap.
    subtraction_expected = (6.0, (7.5, 5.25, 4.5),
                            ((14.0, 6.0, 5.625), (6.0, 6.75, 3.9375),
                             (5.625, 3.9375, 4.8125)))
    subtraction_moments = []
    subtraction_digests = []
    for repetition in (1, 2):
        output = work / f"subtraction-{repetition}.obj"
        subprocess.run((str(executable), "subtraction", str(left), str(right), str(output)),
                       check=True, capture_output=True, text=True, timeout=120)
        moments = integrate_geometric_moments(read_obj(output))
        observed = (moments.volume, moments.first_moment, moments.second_moment)
        subtraction_moments.append(observed)
        subtraction_digests.append(hashlib.sha256(output.read_bytes()).hexdigest())
        for actual, target in zip((observed[0], *observed[1], *(v for row in observed[2] for v in row)),
                                  (subtraction_expected[0], *subtraction_expected[1],
                                   *(v for row in subtraction_expected[2] for v in row))):
            assert abs(actual - target) <= 1e-9 * max(1.0, abs(target)), (actual, target)
        print(f"subtraction_run={repetition} volume={moments.volume:.17g} "
              f"first_moment={moments.first_moment} sha256={subtraction_digests[-1]}")
    assert subtraction_moments[0] == subtraction_moments[1]
    assert subtraction_digests[0] == subtraction_digests[1], subtraction_digests
    print("SYNTHETIC_REGRESSION BOOLEAN_COMPONENT_ONLY PASS "
          "repeat_mass_properties_equal=True raw_output_hashes_equal=True")


if __name__ == "__main__":
    main()
