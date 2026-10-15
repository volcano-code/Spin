"""Parse compile-time-only stage timing events without granting scientific status."""

import re


_EVENT = re.compile(
    r"^SPIN2014_DIAG_TIMING phase=([a-z][a-z0-9_.]*) "
    r"pass=(-?\d+) elapsed_ns=(\d+)$"
)
_PREFIX = "SPIN2014_DIAG_TIMING "


def parse_timing_log(text):
    phases = {}
    event_count = 0
    for line in text.splitlines():
        if not line.startswith(_PREFIX):
            continue
        match = _EVENT.fullmatch(line)
        if match is None:
            raise ValueError("MALFORMED_DIAGNOSTIC_TIMER: " + line)
        phase, pass_index, elapsed_ns = match.groups()
        elapsed_ns = int(elapsed_ns)
        item = phases.setdefault(phase, {"count": 0, "elapsed_ns": 0, "events": []})
        item["count"] += 1
        item["elapsed_ns"] += elapsed_ns
        item["events"].append({"pass": int(pass_index), "elapsed_ns": elapsed_ns})
        event_count += 1
    for item in phases.values():
        item["events"].sort(key=lambda event: event["pass"])
    return {"event_count": event_count, "phases": dict(sorted(phases.items()))}
