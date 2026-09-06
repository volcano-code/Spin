"""Deliberately bounded GNU C/C++ command and make-dependency grammar.

Only one ordinary .c/.cc/.cpp/.cxx/.C source, -c and explicit -o are
accepted. Preprocessed sources, assembler, response files, plugins, specs,
PCH/modules, forwarding options and unrecognized flags fail closed. This is
not a general compiler command parser. Extend admission with behavioral tests.
"""

from dataclasses import dataclass
from pathlib import Path
import re


class ObservationError(ValueError):
    """The invocation cannot be observed under the supported contract."""


@dataclass(frozen=True)
class CompileCommand:
    source: str
    output: str
    depfile: str | None
    probe_argv: list[str]


def parse_command(argv: list[str]) -> CompileCommand:
    if not argv or any(not isinstance(a, str) or not a or "\x00" in a or "\n" in a or "\r" in a for a in argv):
        raise ObservationError("expected nonempty compiler argv without NUL/newline")
    if any(a.startswith("@") for a in argv):
        raise ObservationError("response files are unsupported")
    source, output, depfile, language = [], [], [], []
    probe = [argv[0]]
    compile_count = 0
    dependencies = False
    values = {"-I", "-isystem", "-iquote", "-idirafter", "-D", "-U", "-include", "-imacros", "-x", "-isysroot", "--sysroot"}
    plain = {"-pipe", "-pthread", "-ansi", "-pedantic", "-pedantic-errors", "-nostdinc", "-nostdinc++",
             "-fPIC", "-fpic", "-fPIE", "-fpie", "-fno-PIE", "-fno-pie", "-fno-PIC", "-fno-pic",
             "-frounding-math", "-fno-rounding-math", "-ffast-math", "-fno-fast-math",
             "-fexceptions", "-fno-exceptions", "-frtti", "-fno-rtti", "-fno-omit-frame-pointer",
             "-fomit-frame-pointer", "-ffunction-sections", "-fdata-sections", "-fvisibility=hidden",
             "-fvisibility=default", "-fvisibility-inlines-hidden", "-fwrapv", "-fno-strict-aliasing",
             "-fstrict-aliasing", "-fno-common", "-fcommon", "-fstack-protector", "-fstack-protector-strong",
             "-m32", "-m64", "-msse", "-msse2", "-msse3", "-mssse3", "-msse4.1", "-msse4.2",
             "-mavx", "-mavx2", "-mfma"}
    i = 1
    while i < len(argv):
        a = argv[i]
        if a in {"-o", "-MF", "-MT", "-MQ"} | values:
            if i + 1 >= len(argv) or argv[i + 1].startswith("-"):
                raise ObservationError("missing or ambiguous argument for " + a)
            value = argv[i + 1]
            i += 2
            if a == "-o":
                output.append(value)
            elif a == "-MF":
                depfile.append(value)
            elif a not in {"-MT", "-MQ"}:
                probe.extend([a, value])
                if a == "-x":
                    language.append(value)
            continue
        if a == "-c":
            compile_count += 1
        elif a in {"-MD", "-MMD"}:
            dependencies = True
        elif a == "-MP":
            pass
        elif a.startswith(("-MF", "-MT", "-MQ")) and len(a) > 3:
            if a.startswith("-MF"):
                depfile.append(a[3:])
        elif a.startswith("-o") and len(a) > 2:
            output.append(a[2:])
        elif a.startswith(("-I", "-D", "-U")) and len(a) > 2:
            probe.append(a)
        elif a.startswith(("-isystem", "-iquote", "-idirafter")) and a not in values:
            probe.append(a)
        elif a.startswith("--sysroot=") and len(a) > len("--sysroot="):
            probe.append(a)
        elif a in plain or re.fullmatch(r"-O(?:[0-3sgz]|fast)?|-g(?:[0-3]|gdb[0-3]?|dwarf-[2-5])?|-W(?:no-)?[A-Za-z0-9_-]+(?:=[A-Za-z0-9_-]+)?|-std=(?:gnu|c)(?:\+\+)?[0-9a-z]+", a):
            probe.append(a)
        elif not a.startswith("-") and Path(a).suffix in {".c", ".cc", ".cpp", ".cxx", ".C"}:
            source.append(a)
            probe.append(a)
        else:
            raise ObservationError("unsupported compiler argument: " + a)
        i += 1
    if compile_count != 1 or len(source) != 1 or len(output) != 1:
        raise ObservationError("require exactly one -c, C/C++ source and explicit -o output")
    if len(depfile) > 1 or (depfile and not dependencies) or (language and language not in (["c"], ["c++"])):
        raise ObservationError("ambiguous dependency output or unsupported language selection")
    if output[0] == "-" or (depfile and depfile[0] == "-"):
        raise ObservationError("stdout outputs are unsupported")
    return CompileCommand(source[0], output[0], depfile[0] if depfile else str(Path(output[0]).with_suffix(".d")) if dependencies else None,
                          probe + ["-M", "-fpch-preprocess", "-MT", "observer-target"])


def parse_dependencies(text: str) -> list[str]:
    """Read one GNU -M rule, including continuations, escaped space/# and $$.

    GNU quotes backslashes only when they precede whitespace; literal
    backslashes preceding normal path characters must be retained. Newlines,
    tabs, colons and ambiguous make syntax in filenames are not supported.
    """
    text = text.replace("\\\n", "")
    prefix = "observer-target:"
    if not text.startswith(prefix):
        raise ObservationError("unexpected dependency rule target")
    text = text[len(prefix):].strip()
    tokens, token = [], []
    i = 0
    while i < len(text):
        c = text[i]
        if c == "\\":
            start = i
            while i < len(text) and text[i] == "\\":
                i += 1
            count = i - start
            if i < len(text) and text[i] in " \t":
                token.extend("\\" * (count // 2))
                if count % 2:
                    token.append(text[i])
                    i += 1
            elif i < len(text) and text[i] == "#":
                token.extend("\\" * (count - 1))
                token.append("#")
                i += 1
            else:
                token.extend("\\" * count)
            continue
        if c.isspace():
            if c not in " \t":
                raise ObservationError("multiple dependency rules or newline paths unsupported")
            if token:
                tokens.append("".join(token))
                token = []
        elif c == "$":
            if i + 1 >= len(text) or text[i + 1] != "$":
                raise ObservationError("unescaped make expansion in dependencies")
            token.append("$")
            i += 1
        elif c in "#:|;":
            raise ObservationError("ambiguous make syntax in dependency path")
        else:
            token.append(c)
        i += 1
    if token:
        tokens.append("".join(token))
    if not tokens:
        raise ObservationError("empty dependency closure")
    return sorted(set(tokens))
