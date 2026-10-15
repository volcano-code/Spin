from stage_timing_audit import parse_timing_log


def test_stage_timer_parser_keeps_phase_and_pass_identity():
    log = "\n".join((
        "worker diagnostic noise",
        "SPIN2014_DIAG_TIMING phase=checkpoint_decode_history pass=18 elapsed_ns=120",
        "SPIN2014_DIAG_TIMING phase=checkpoint_decode_history pass=19 elapsed_ns=280",
        "SPIN2014_DIAG_TIMING phase=adaptive_split pass=20 elapsed_ns=900",
    ))

    parsed = parse_timing_log(log)

    assert parsed["event_count"] == 3
    assert parsed["phases"]["checkpoint_decode_history"] == {
        "count": 2,
        "elapsed_ns": 400,
        "events": [
            {"pass": 18, "elapsed_ns": 120},
            {"pass": 19, "elapsed_ns": 280},
        ],
    }
    assert parsed["phases"]["adaptive_split"]["events"] == [
        {"pass": 20, "elapsed_ns": 900},
    ]


def test_stage_timer_parser_rejects_malformed_timer_events():
    try:
        parse_timing_log("SPIN2014_DIAG_TIMING phase=adaptive_split elapsed_ns=not-a-number")
    except ValueError as exc:
        assert "MALFORMED_DIAGNOSTIC_TIMER" in str(exc)
    else:
        raise AssertionError("malformed diagnostic timing record was accepted")
