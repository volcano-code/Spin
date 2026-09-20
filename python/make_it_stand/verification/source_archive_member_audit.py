"""Bounded, static GNU-ar member-byte check; never a build execution attestation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys


class ArchiveMemberAuditError(ValueError):
    """Evidence is malformed, ambiguous, unsupported, or no longer matches disk."""


_SHA_LEN = 64
_AR_MAGIC = b"!<arch>\n"
_THIN_AR_MAGIC = b"!<thin>\n"


def _fail(message: str) -> None:
    raise ArchiveMemberAuditError(message)


def _mapping(value, label):
    if type(value) is not dict:
        _fail(label + " must be an object")
    return value


def _list(value, label):
    if type(value) is not list:
        _fail(label + " must be an array")
    return value


def _string(value, label):
    if type(value) is not str or not value:
        _fail(label + " must be a nonempty string")
    return value


def _sha(value, label):
    value = _string(value, label)
    if len(value) != _SHA_LEN or any(c not in "0123456789abcdef" for c in value):
        _fail(label + " must be a lowercase SHA-256")
    return value


def _size(value, label):
    if type(value) is not int or value < 0:
        _fail(label + " must be a nonnegative integer")
    return value


def _load_json(path, label):
    try:
        with Path(path).open("r", encoding="utf-8") as stream:
            return _mapping(json.load(stream), label)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        _fail("cannot read " + label + ": " + str(error))


def _stable_hash(path, label):
    path = Path(path)
    try:
        before = path.stat()
        if not stat.S_ISREG(before.st_mode):
            _fail(label + " is not a regular file")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            opened = os.fstat(stream.fileno())
            if not stat.S_ISREG(opened.st_mode):
                _fail(label + " is not a regular file")
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
            closed = os.fstat(stream.fileno())
        after = path.stat()
    except OSError as error:
        _fail("cannot hash " + label + ": " + str(error))
    identity = lambda item: (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
    if identity(before) != identity(opened) or identity(opened) != identity(closed) or identity(closed) != identity(after):
        _fail(label + " changed while hashing")
    return {"path": str(path), "resolved_path": str(path.resolve(strict=True)),
            "size": after.st_size, "sha256": digest.hexdigest()}


def _check_file_record(record, label, *, require_resolved=False):
    record = _mapping(record, label)
    path = Path(_string(record.get("path"), label + ".path"))
    observed = _stable_hash(path, label)
    if observed["size"] != _size(record.get("size"), label + ".size") or observed["sha256"] != _sha(record.get("sha256"), label + ".sha256"):
        _fail(label + " does not match its recorded bytes")
    if require_resolved and observed["resolved_path"] != _string(record.get("resolved_path"), label + ".resolved_path"):
        _fail(label + " resolved path differs from inventory")
    return observed


def _same_file_bytes(left, right):
    return all(left[key] == right[key] for key in ("resolved_path", "size", "sha256"))


def _safe_relative_object(value, label):
    value = _string(value, label)
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in ("", ".", "..") for part in path.parts):
        _fail(label + " is not a safe relative object path")
    if any("\\" in part or "\x00" in part for part in path.parts) or path.suffix != ".o":
        _fail(label + " is not a supported object argument")
    return value


def _safe_member_name(value, label):
    value = _string(value, label)
    if (value in (".", "..") or "/" in value or "\\" in value or "\x00" in value or
            value != Path(value).name or not value.endswith(".o")):
        _fail(label + " is unsafe or unsupported")
    return value


def _parse_static_archive_command(command):
    """Accept only the two observed CMake/Ninja static-ar command spellings."""
    command = _string(command, "generated archive command")
    if any(mark in command for mark in ("\n", "\r", "`", "$", "|", ";", "<", ">", "'", '"')):
        _fail("generated archive command uses unsupported shell syntax")
    tokens = command.split()
    prefix = [":", "&&", "/usr/bin/cmake", "-E", "rm", "-f"]
    if tokens[:6] != prefix or len(tokens) < 12:
        _fail("generated archive command is not the observed GNU/CMake form")
    removed = tokens[6]
    if tokens[7:10] != ["&&", "/usr/bin/ar", "qc"]:
        _fail("generated archive command is not the observed GNU/CMake form")
    start = 10
    if len(tokens) < start + 5:
        _fail("generated archive command has no member arguments")
    archive = tokens[start]
    if PurePosixPath(archive).is_absolute() or archive in (".", "..") or ".." in PurePosixPath(archive).parts:
        _fail("generated archive output is unsafe")
    if removed is not None and removed != archive:
        _fail("archive removal and archive output differ")
    try:
        stop = tokens.index("&&", start + 1)
    except ValueError:
        _fail("generated archive command has no bounded ar terminator")
    objects = [_safe_relative_object(item, "archive object") for item in tokens[start + 1:stop]]
    if not objects or tokens[stop:stop + 5] != ["&&", "/usr/bin/ranlib", archive, "&&", ":"] or len(tokens) != stop + 5:
        _fail("generated archive command has an unsupported tail")
    return archive, objects


def _parse_gnu_ar_members(archive):
    """Read ordinary GNU ar directly, avoiding extraction and any evidence writes."""
    archive = Path(archive)
    with archive.open("rb") as stream:
        magic = stream.read(8)
        if magic == _THIN_AR_MAGIC:
            _fail("thin archives are unsupported")
        if magic != _AR_MAGIC:
            _fail("archive is not ordinary GNU ar")
        end = archive.stat().st_size
        string_table = None
        result = []
        seen = set()
        while stream.tell() < end:
            header = stream.read(60)
            if len(header) != 60 or header[58:60] != b"`\n":
                _fail("malformed GNU ar member header")
            raw_size = header[48:58].strip()
            if not raw_size or not raw_size.isdigit():
                _fail("malformed GNU ar member size")
            size = int(raw_size)
            payload_offset = stream.tell()
            if size > end - payload_offset:
                _fail("GNU ar member exceeds archive bounds")
            raw_name = header[:16].rstrip(b" ")
            if raw_name == b"//":
                if string_table is not None or size > 16 * 1024 * 1024:
                    _fail("unsupported GNU ar long-name table")
                string_table = stream.read(size)
            elif raw_name == b"/":
                stream.seek(size, os.SEEK_CUR)  # GNU symbol table, never a material member.
            else:
                if raw_name.startswith(b"#1/"):
                    _fail("BSD extended ar names are unsupported")
                if raw_name.startswith(b"/"):
                    if string_table is None or not raw_name[1:].isdigit():
                        _fail("unsupported GNU ar member name")
                    offset = int(raw_name[1:])
                    if offset >= len(string_table):
                        _fail("GNU ar name offset is out of range")
                    marker = string_table.find(b"/\n", offset)
                    if marker < 0:
                        _fail("unterminated GNU ar long member name")
                    name_bytes = string_table[offset:marker]
                else:
                    if not raw_name.endswith(b"/"):
                        _fail("unsupported GNU ar short member name")
                    name_bytes = raw_name[:-1]
                try:
                    name = name_bytes.decode("ascii", "strict")
                except UnicodeDecodeError:
                    _fail("GNU ar member name is not ASCII")
                _safe_member_name(name, "archive member")
                if name in seen:
                    _fail("duplicate archive member name: " + name)
                seen.add(name)
                result.append({"name": name, "offset": payload_offset, "size": size})
                stream.seek(size, os.SEEK_CUR)
            if stream.tell() & 1:
                if stream.read(1) != b"\n":
                    _fail("malformed GNU ar member padding")
        if stream.tell() != end:
            _fail("GNU ar parser did not consume the archive")
    return result


def _hash_member(archive, offset, size):
    digest = hashlib.sha256()
    with Path(archive).open("rb") as stream:
        stream.seek(offset)
        remaining = size
        while remaining:
            block = stream.read(min(1024 * 1024, remaining))
            if not block:
                _fail("archive member became truncated during inspection")
            digest.update(block)
            remaining -= len(block)
    return digest.hexdigest()


def _ar_identity_and_listing(archive, expected_names):
    ar = Path("/usr/bin/ar")
    observed = _stable_hash(ar, "/usr/bin/ar")
    try:
        completed = subprocess.run([str(ar), "t", str(archive)], check=False, capture_output=True,
                                   timeout=30, text=True, encoding="ascii", errors="strict")
    except (OSError, subprocess.SubprocessError, UnicodeError) as error:
        _fail("/usr/bin/ar listing failed: " + str(error))
    if completed.returncode != 0 or completed.stderr:
        _fail("/usr/bin/ar could not list the archive")
    listed = completed.stdout.splitlines()
    if listed != expected_names:
        _fail("/usr/bin/ar listing disagrees with direct archive parsing")
    return observed


def _output_argument(argv):
    argv = _list(argv, "compiler receipt argv")
    places = [index for index, item in enumerate(argv) if item == "-o"]
    if len(places) != 1 or places[0] + 1 >= len(argv):
        _fail("compiler receipt does not have one explicit -o output")
    return _string(argv[places[0] + 1], "compiler receipt output")


def _bind_compiler_receipt(inventory_receipt, expected_logical):
    inventory_receipt = _mapping(inventory_receipt, "compiler receipt inventory entry")
    if inventory_receipt.get("classification") != "target_compilation" or inventory_receipt.get("status") != "OBSERVED_SUCCESS":
        _fail("archive object does not have an observed successful target compiler receipt")
    receipt_record = _mapping(inventory_receipt.get("receipt"), "compiler receipt record")
    _check_file_record(receipt_record, "compiler receipt")
    receipt_path = Path(_string(receipt_record.get("path"), "compiler receipt.path"))
    raw = _load_json(receipt_path, "compiler receipt")
    if raw.get("schema") != "make-it-stand-compiler-observation-v1" or raw.get("status") != "OBSERVED_SUCCESS" or raw.get("compiler_exit_status") != 0:
        _fail("raw compiler receipt is not an observed successful compiler receipt")
    if _output_argument(raw.get("argv")) != expected_logical:
        _fail("compiler receipt output argument does not match generated archive command")
    cwd = Path(_string(raw.get("cwd"), "compiler receipt cwd"))
    recorded = _mapping(inventory_receipt.get("recorded_object"), "recorded object")
    actual = _mapping(inventory_receipt.get("actual_object"), "actual object")
    raw_object = _mapping(raw.get("object"), "raw compiler receipt object")
    for label, object_record in (("recorded object", recorded), ("raw compiler receipt object", raw_object)):
        if object_record.get("logical_path") != expected_logical:
            _fail(label + " logical output differs from generated archive command")
        _size(object_record.get("size"), label + ".size")
        _sha(object_record.get("sha256"), label + ".sha256")
        _string(object_record.get("resolved_path"), label + ".resolved_path")
    if str(cwd / expected_logical) != raw_object["resolved_path"]:
        _fail("raw compiler receipt cwd/output identity is inconsistent")
    actual_path = Path(_string(actual.get("path"), "actual object.path"))
    observed = _stable_hash(actual_path, "actual compiler object")
    for label, object_record in (("recorded object", recorded), ("raw compiler receipt object", raw_object), ("actual object", actual)):
        if object_record.get("size") != observed["size"] or object_record.get("sha256") != observed["sha256"]:
            _fail(label + " bytes do not match the actual compiler object")
    if str(actual_path) != recorded["resolved_path"] or recorded["resolved_path"] != raw_object["resolved_path"]:
        _fail("compiler object paths disagree")
    return {"name": _safe_member_name(PurePosixPath(expected_logical).name, "expected archive member"),
            "path": str(actual_path), "size": observed["size"], "sha256": observed["sha256"]}


def audit_archive(inventory_path, ninja_evidence_path, archive_path):
    """Bind one captured static archive to observed compiler-object bytes.

    The return value deliberately leaves execution, full build provenance, runtime,
    solver, paper, and physics acceptance out of scope.
    """
    inventory = _load_json(inventory_path, "build inventory")
    evidence = _load_json(ninja_evidence_path, "Ninja command evidence")
    archive_path = Path(archive_path)
    observed_archive = _stable_hash(archive_path, "requested archive")
    archives = _list(inventory.get("archives"), "build inventory archives")
    archive_records = [record for record in archives if type(record) is dict and
                       str(Path(record.get("path", "")).resolve(strict=False)) == observed_archive["resolved_path"]]
    if len(archive_records) != 1:
        _fail("requested archive has no unique inventory record")
    recorded_archive = _check_file_record(archive_records[0], "inventory archive")
    if recorded_archive["resolved_path"] != observed_archive["resolved_path"]:
        _fail("inventory archive path differs from requested archive")
    source_record = _mapping(evidence.get("command_source"), "Ninja command source")
    _check_file_record({"path": source_record.get("resolved_path"), "size": source_record.get("size"),
                        "sha256": source_record.get("sha256")}, "Ninja command source")
    commands = _list(evidence.get("generated_archive_and_link_commands"), "generated archive/link commands")
    selected = []
    for entry in commands:
        entry = _mapping(entry, "generated archive/link command")
        if entry.get("kind") != "static_archive":
            continue
        archive_argument, objects = _parse_static_archive_command(entry.get("command"))
        if observed_archive["resolved_path"].endswith("/" + archive_argument):
            selected.append((archive_argument, objects))
    if len(selected) != 1:
        _fail("requested archive has no unique captured static archive command")
    _, object_arguments = selected[0]
    receipts = _list(inventory.get("compiler_receipts"), "build inventory compiler receipts")
    expected = []
    for argument in object_arguments:
        matches = [entry for entry in receipts if type(entry) is dict and
                   type(entry.get("recorded_object")) is dict and entry["recorded_object"].get("logical_path") == argument]
        if len(matches) != 1:
            _fail("archive object does not have exactly one inventory compiler receipt: " + argument)
        expected.append(_bind_compiler_receipt(matches[0], argument))
    expected_names = [entry["name"] for entry in expected]
    if len(set(expected_names)) != len(expected_names):
        _fail("generated archive command has duplicate object basenames")
    members = _parse_gnu_ar_members(archive_path)
    member_names = [entry["name"] for entry in members]
    ar_identity = _ar_identity_and_listing(archive_path, member_names)
    if member_names != expected_names:
        _fail("archive member names are missing, extra, duplicated, or reordered")
    member_bytes = 0
    bound_members = []
    for member, object_record in zip(members, expected):
        digest = _hash_member(archive_path, member["offset"], member["size"])
        if member["size"] != object_record["size"] or digest != object_record["sha256"]:
            _fail("archive member bytes differ from compiler object: " + member["name"])
        member_bytes += member["size"]
        bound_members.append({"name": member["name"], "size": member["size"], "sha256": digest,
                              "compiler_object": object_record["path"]})
    final_archive = _stable_hash(archive_path, "archive after member inspection")
    if not _same_file_bytes(final_archive, observed_archive) or not _same_file_bytes(final_archive, recorded_archive):
        _fail("archive changed during verification or no longer matches inventory")
    return {"schema": "make-it-stand-source-archive-member-audit-v1",
            "status": "OBSERVED_ARCHIVE_MEMBER_BYTES_MATCH",
            "archive": final_archive, "ar_executable": ar_identity,
            "member_count": len(bound_members), "member_bytes": member_bytes,
            "members": bound_members,
            "archive_command_execution": "NOT_VERIFIED",
            "build_source_attestation": "NOT_VERIFIED",
            "full_build_runtime_attestation": "NOT_VERIFIED",
            "solver_acceptance": "NOT_VERIFIED", "paper_acceptance": "NOT_VERIFIED",
            "physics_acceptance": "NOT_VERIFIED",
            "limitations": ["Ninja generated-command evidence is not an archive execution trace",
                            "only ordinary GNU ar member bytes are bound to recorded compiler object bytes",
                            "linker implicit inputs, compiler/assembler plugins, runtime bytes, source, solver, paper and physics remain unverified"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--ninja-evidence", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        result = audit_archive(args.inventory, args.ninja_evidence, args.archive)
        with Path(args.output).open("x", encoding="utf-8") as stream:
            json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write("\n")
    except (ArchiveMemberAuditError, OSError) as error:
        print("archive member audit failed: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
