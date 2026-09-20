"""Real GNU-ar fixtures for the bounded archive-member provenance check."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pytest


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(path: Path) -> dict:
    return {"path": str(path), "size": path.stat().st_size, "sha256": _sha(path)}


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _fixture(tmp_path: Path, names=("first.o", "second.o")):
    """Make an archive from raw bytes, plus independently hand-built receipts."""
    ar = shutil.which("ar")
    assert ar, "real GNU ar is required for this behavioral fixture"
    build = tmp_path / "build"
    objects = build / "objects"
    objects.mkdir(parents=True)
    members = []
    receipts = []
    for index, name in enumerate(names):
        obj = objects / name
        obj.parent.mkdir(parents=True, exist_ok=True)
        obj.write_bytes((f"object-{index}-bytes\n").encode())
        logical = str(obj.relative_to(build))
        receipt = build / "receipts" / f"{index}.json"
        receipt.parent.mkdir(exist_ok=True)
        receipt_value = {
            "schema": "make-it-stand-compiler-observation-v1",
            "status": "OBSERVED_SUCCESS",
            "compiler_exit_status": 0,
            "cwd": str(build),
            "argv": ["/usr/bin/cc", "-c", f"source-{index}.c", "-o", logical],
            "object": {"logical_path": logical, "resolved_path": str(obj),
                       "size": obj.stat().st_size, "sha256": _sha(obj)},
        }
        _write_json(receipt, receipt_value)
        receipts.append({
            "classification": "target_compilation", "status": "OBSERVED_SUCCESS",
            "receipt": _record(receipt),
            "recorded_object": {"logical_path": logical, "resolved_path": str(obj),
                                "size": obj.stat().st_size, "sha256": _sha(obj)},
            "actual_object": _record(obj),
        })
        members.append(obj)
    archive = build / "libfixture.a"
    subprocess.run([ar, "qc", str(archive), *map(str, members)], check=True)
    subprocess.run([ar, "s", str(archive)], check=True)
    ninja_output = tmp_path / "ninja-commands.stdout"
    command_members = " ".join(str(p.relative_to(build)) for p in members)
    ninja_output.write_text(
        f": && /usr/bin/cmake -E rm -f libfixture.a && /usr/bin/ar qc libfixture.a "
        f"{command_members} && /usr/bin/ranlib libfixture.a && :\n"
    )
    inventory = tmp_path / "inventory.json"
    _write_json(inventory, {"archives": [_record(archive)], "compiler_receipts": receipts})
    evidence = tmp_path / "ninja-link-evidence.json"
    _write_json(evidence, {
        "command_source": {"logical_path": str(ninja_output), "resolved_path": str(ninja_output),
                           "size": ninja_output.stat().st_size, "sha256": _sha(ninja_output)},
        "generated_archive_and_link_commands": [{"kind": "static_archive", "line_number": 1,
                                                   "command": ninja_output.read_text().strip()}],
    })
    return {"archive": archive, "inventory": inventory, "evidence": evidence,
            "objects": members, "ninja_output": ninja_output, "ar": ar}


def _audit(case):
    from make_it_stand.verification.source_archive_member_audit import audit_archive
    return audit_archive(case["inventory"], case["evidence"], case["archive"])


def _failure(case):
    from make_it_stand.verification.source_archive_member_audit import ArchiveMemberAuditError
    with pytest.raises(ArchiveMemberAuditError):
        _audit(case)


def test_real_gnu_archive_binds_exact_ordered_member_bytes(tmp_path):
    case = _fixture(tmp_path)
    result = _audit(case)
    assert result["status"] == "OBSERVED_ARCHIVE_MEMBER_BYTES_MATCH"
    assert result["archive_command_execution"] == "NOT_VERIFIED"
    assert result["full_build_runtime_attestation"] == "NOT_VERIFIED"
    assert result["member_count"] == 2
    assert result["member_bytes"] == sum(p.stat().st_size for p in case["objects"])
    ar_identity = result["ar_executable"]
    assert ar_identity["resolved_path"].startswith("/usr/bin/")
    assert "ar" in Path(ar_identity["resolved_path"]).name


def test_same_member_names_with_altered_payload_fails_closed(tmp_path):
    case = _fixture(tmp_path)
    altered = tmp_path / "altered" / "first.o"
    altered.parent.mkdir()
    altered.write_bytes(b"changed object payload\n")
    subprocess.run([case["ar"], "r", str(case["archive"]), str(altered)], check=True)
    inventory = json.loads(case["inventory"].read_text())
    inventory["archives"] = [_record(case["archive"])]
    _write_json(case["inventory"], inventory)
    _failure(case)


def test_thin_archive_is_rejected(tmp_path):
    case = _fixture(tmp_path)
    case["archive"].unlink()
    subprocess.run([case["ar"], "rcT", str(case["archive"]), *map(str, case["objects"])], check=True)
    inventory = json.loads(case["inventory"].read_text())
    inventory["archives"] = [_record(case["archive"])]
    _write_json(case["inventory"], inventory)
    _failure(case)


def test_duplicate_member_or_expected_basename_is_rejected(tmp_path):
    case = _fixture(tmp_path, names=("left/same.o", "right/same.o"))
    _failure(case)


@pytest.mark.parametrize("mode", ("missing", "extra", "reordered", "unsafe"))
def test_wrong_or_unsafe_member_sequence_is_rejected(tmp_path, mode):
    case = _fixture(tmp_path)
    if mode == "missing":
        subprocess.run([case["ar"], "d", str(case["archive"]), "second.o"], check=True)
    elif mode == "extra":
        extra = tmp_path / "extra.o"; extra.write_bytes(b"extra\n")
        subprocess.run([case["ar"], "r", str(case["archive"]), str(extra)], check=True)
    elif mode == "reordered":
        subprocess.run([case["ar"], "d", str(case["archive"]), "first.o"], check=True)
        subprocess.run([case["ar"], "q", str(case["archive"]), str(case["objects"][0])], check=True)
    else:
        data = case["archive"].read_bytes()
        header = bytearray(data[8:68]); header[:16] = b"../unsafe.o/    "
        case["archive"].write_bytes(data[:8] + bytes(header) + data[68:])
    inventory = json.loads(case["inventory"].read_text())
    inventory["archives"] = [_record(case["archive"])]
    _write_json(case["inventory"], inventory)
    _failure(case)


@pytest.mark.parametrize("mode", ("receipt", "object"))
def test_malformed_receipt_or_object_identity_mismatch_is_rejected(tmp_path, mode):
    case = _fixture(tmp_path)
    inventory = json.loads(case["inventory"].read_text())
    if mode == "receipt":
        receipt = Path(inventory["compiler_receipts"][0]["receipt"]["path"])
        receipt.write_text("{ malformed\n")
        inventory["compiler_receipts"][0]["receipt"] = _record(receipt)
    else:
        inventory["compiler_receipts"][0]["recorded_object"]["logical_path"] = "objects/not-first.o"
    _write_json(case["inventory"], inventory)
    _failure(case)


def test_digest_mismatch_and_archive_mutation_are_rejected(tmp_path):
    case = _fixture(tmp_path)
    inventory = json.loads(case["inventory"].read_text())
    inventory["archives"][0]["sha256"] = "0" * 64
    _write_json(case["inventory"], inventory)
    _failure(case)


def test_archive_mutation_after_member_inspection_is_rejected(tmp_path, monkeypatch):
    case = _fixture(tmp_path)
    from make_it_stand.verification import source_archive_member_audit as api
    original = api._hash_member
    calls = 0

    def mutate_after_last_member(*args):
        nonlocal calls
        result = original(*args)
        calls += 1
        if calls == len(case["objects"]):
            data = case["archive"].read_bytes()
            case["archive"].write_bytes(data[:-1] + bytes([data[-1] ^ 1]))
        return result

    monkeypatch.setattr(api, "_hash_member", mutate_after_last_member)
    with pytest.raises(api.ArchiveMemberAuditError, match="changed during verification"):
        api.audit_archive(case["inventory"], case["evidence"], case["archive"])


def test_structured_archive_command_must_match_its_hashed_raw_ninja_line(tmp_path):
    case = _fixture(tmp_path)
    case["ninja_output"].write_text("unrelated generated command\n")
    evidence = json.loads(case["evidence"].read_text())
    evidence["command_source"] = {
        "logical_path": str(case["ninja_output"]), "resolved_path": str(case["ninja_output"]),
        "size": case["ninja_output"].stat().st_size, "sha256": _sha(case["ninja_output"]),
    }
    _write_json(case["evidence"], evidence)
    _failure(case)


def test_carriage_return_cannot_invent_a_raw_ninja_command_line(tmp_path):
    case = _fixture(tmp_path)
    evidence = json.loads(case["evidence"].read_text())
    command = evidence["generated_archive_and_link_commands"][0]
    case["ninja_output"].write_bytes(b"unrelated\r" + command["command"].encode() + b"\n")
    command["line_number"] = 2
    evidence["command_source"] = {
        "logical_path": str(case["ninja_output"]), "resolved_path": str(case["ninja_output"]),
        "size": case["ninja_output"].stat().st_size, "sha256": _sha(case["ninja_output"]),
    }
    _write_json(case["evidence"], evidence)
    _failure(case)


@pytest.mark.parametrize("mode", ("invalid_line", "duplicate_line", "ambiguous_target"))
def test_invalid_duplicate_or_ambiguous_static_command_records_fail_closed(tmp_path, mode):
    case = _fixture(tmp_path)
    evidence = json.loads(case["evidence"].read_text())
    command = evidence["generated_archive_and_link_commands"][0]
    if mode == "invalid_line":
        command["line_number"] = 0
    elif mode == "duplicate_line":
        evidence["generated_archive_and_link_commands"].append(dict(command))
    else:
        case["ninja_output"].write_text(case["ninja_output"].read_text() * 2)
        second = dict(command)
        second["line_number"] = 2
        evidence["generated_archive_and_link_commands"].append(second)
        evidence["command_source"] = {
            "logical_path": str(case["ninja_output"]), "resolved_path": str(case["ninja_output"]),
            "size": case["ninja_output"].stat().st_size, "sha256": _sha(case["ninja_output"]),
        }
    _write_json(case["evidence"], evidence)
    _failure(case)
