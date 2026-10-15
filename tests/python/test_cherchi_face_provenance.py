from __future__ import annotations

from pathlib import Path

from make_it_stand.verification.cherchi_face_provenance import recover_face_sources


def write_triangle(path: Path) -> None:
    path.write_text(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n", encoding="ascii"
    )


def exact_rows() -> list[list[str]]:
    return [
        ["0", "MATCHED", "0/1", "0/1", "0/1", "0", "0", "0", "0", "0", "0"],
        ["1", "MATCHED", "0b1/0b1", "0/1", "0/1", "1", "0", "0", "1", "0", "0"],
        ["2", "MATCHED", "0/1", "0b1/0b1", "0/1", "0", "1", "0", "0", "1", "0"],
    ]


def test_coincident_input_facets_are_reported_as_ambiguous(tmp_path: Path) -> None:
    first, second, output = (tmp_path / "first.obj", tmp_path / "second.obj",
                             tmp_path / "output.obj")
    write_triangle(first)
    write_triangle(second)
    write_triangle(output)

    records = recover_face_sources(
        (first, second), output, exact_rows(), [["0", "0x1"]]
    )

    assert records == [{
        "face_index": 0,
        "input_object_mask": 0x1,
        "source_triangles": [
            {"input_object": 0, "input_triangle": 0},
            {"input_object": 1, "input_triangle": 0},
        ],
        "status": "AMBIGUOUS",
    }]


def test_face_outside_all_input_facets_fails_closed_as_unmatched(tmp_path: Path) -> None:
    source, output = tmp_path / "source.obj", tmp_path / "output.obj"
    write_triangle(source)
    output.write_text(
        "v 0 0 1\nv 1 0 1\nv 0 1 1\nf 1 2 3\n", encoding="ascii"
    )
    shifted_rows = [
        ["0", "MATCHED", "0/1", "0/1", "0b1/0b1", "0", "0", "1", "0", "0", "1"],
        ["1", "MATCHED", "0b1/0b1", "0/1", "0b1/0b1", "1", "0", "1", "1", "0", "1"],
        ["2", "MATCHED", "0/1", "0b1/0b1", "0b1/0b1", "0", "1", "1", "0", "1", "1"],
    ]

    records = recover_face_sources((source,), output, shifted_rows,
                                   [["0", "0x1"]])

    assert records[0]["status"] == "UNMATCHED"
    assert records[0]["input_object_mask"] == 0x1
    assert records[0]["source_triangles"] == []
