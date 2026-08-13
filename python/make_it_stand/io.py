"""Deterministic, atomic file helpers for the R9 Python layer."""

from __future__ import annotations

import csv
import io as string_io
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Any, Iterable, Mapping, Sequence


class DataFormatError(ValueError):
    """Raised when JSON input is ambiguous or not standards compliant."""


class UnsafePathError(ValueError):
    """Raised when an output path would traverse a symbolic link."""


def _reject_json_constant(value: str) -> None:
    raise DataFormatError(f"non-finite JSON number is not allowed: {value}")


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DataFormatError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def loads_json_strict(text: str) -> Any:
    """Parse JSON while rejecting duplicate keys and NaN/Infinity."""

    try:
        return json.loads(
            text,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_json_constant,
        )
    except json.JSONDecodeError as error:
        raise DataFormatError(
            f"invalid JSON at line {error.lineno}, column {error.colno}: "
            f"{error.msg}"
        ) from error


def read_json(path: os.PathLike[str] | str) -> Any:
    source = Path(path)
    try:
        text = source.read_text(encoding="utf-8")
    except UnicodeDecodeError as error:
        raise DataFormatError(f"JSON is not valid UTF-8: {source}") from error
    return loads_json_strict(text)


def prepare_output_root(path: os.PathLike[str] | str) -> Path:
    """Create and resolve the caller-selected output root.

    The root itself is an explicit CLI/API argument and may intentionally be a
    symbolic link.  Once selected, its resolved directory becomes the trusted
    anchor; schema-controlled relative paths may not traverse further links.
    """

    selected = Path(path)
    selected.mkdir(parents=True, exist_ok=True)
    root = selected.resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"output root is not a directory: {selected}")
    return root


def _ensure_relative_directory_without_symlinks(root: Path, relative: Path) -> Path:
    current = root
    for component in relative.parts:
        current /= component
        try:
            metadata = os.lstat(current)
        except FileNotFoundError:
            try:
                os.mkdir(current)
            except FileExistsError:
                # Another process may have populated this component between
                # lstat and mkdir.  The lstat below still fails closed on a
                # link or non-directory.
                pass
            metadata = os.lstat(current)
        if stat.S_ISLNK(metadata.st_mode):
            raise UnsafePathError(
                f"refusing symbolic-link output parent below output root: {current}"
            )
        if not stat.S_ISDIR(metadata.st_mode):
            raise NotADirectoryError(
                f"output directory component is not a directory: {current}"
            )
    return current


def prepare_output_path(
    output_directory: os.PathLike[str] | str,
    relative_path: os.PathLike[str] | str,
) -> Path:
    """Return a safe output file path rooted below *output_directory*.

    This is a second line of defence behind configuration validation.  It
    rejects traversal, symbolic-link roots/parents, and a pre-existing link or
    directory at the final filename while retaining nested output folders.
    """

    relative = Path(relative_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise UnsafePathError(
            "output path must be relative and must not contain parent traversal"
        )
    if relative.as_posix() in {"", "."}:
        raise UnsafePathError("output path must name a file")

    root = prepare_output_root(output_directory)
    destination = root / relative
    _ensure_relative_directory_without_symlinks(root, relative.parent)
    try:
        metadata = os.lstat(destination)
    except FileNotFoundError:
        return destination
    if stat.S_ISLNK(metadata.st_mode):
        raise UnsafePathError(f"refusing symbolic-link output file: {destination}")
    if stat.S_ISDIR(metadata.st_mode):
        raise IsADirectoryError(f"output filename is an existing directory: {destination}")
    return destination


def _atomic_replace_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            temporary.write(payload)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def write_text_atomic(
    path: os.PathLike[str] | str,
    text: str,
    *,
    encoding: str = "utf-8",
) -> None:
    if not isinstance(text, str):
        raise TypeError("text must be a string")
    _atomic_replace_bytes(Path(path), text.encode(encoding))


def write_json_atomic(path: os.PathLike[str] | str, value: Any) -> None:
    """Write canonical human-readable JSON without host-dependent metadata."""

    try:
        payload = json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise DataFormatError(f"value cannot be encoded as strict JSON: {error}") from error
    write_text_atomic(path, payload + "\n")


def write_csv_atomic(
    path: os.PathLike[str] | str,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> None:
    """Write a deterministic RFC-4180-style CSV with an explicit header."""

    names = tuple(fieldnames)
    if not names or any(not isinstance(name, str) or not name for name in names):
        raise ValueError("fieldnames must be non-empty strings")
    if len(set(names)) != len(names):
        raise ValueError("fieldnames must be unique")

    buffer = string_io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=names,
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    write_text_atomic(path, buffer.getvalue())


__all__ = [
    "DataFormatError",
    "UnsafePathError",
    "loads_json_strict",
    "prepare_output_root",
    "prepare_output_path",
    "read_json",
    "write_csv_atomic",
    "write_json_atomic",
    "write_text_atomic",
]
