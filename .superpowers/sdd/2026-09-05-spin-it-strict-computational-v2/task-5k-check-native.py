"""Independent original-stencil checks of actual native BBW probe output.

This checker imports only the independent reference, never the project/native
FEM or BBW implementation. Literal rectangular fixtures are synthetic tests.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("bbw_reference", HERE / "task-5k-reference-bbw.py")
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)


def load(path):
    return json.loads(Path(path).read_text(), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def check(native, expected):
    for key in ("cells_per_axis", "rule", "protocol"):
        assert native[key] == expected[key], key
    assert native["scope"] == "NATIVE_SYNTHETIC_BBW_NOT_PAPER_CASE"
    data = oracle.reference(native["cells_per_axis"], native["rule"])
    nodes = np.asarray(native["nodes"], dtype=float)
    assert np.array_equal(nodes, data["nodes"]), "native represented nodes changed"
    # Match actual complete handles by node identity, not assumed column order.
    handles = native["handle_nodes"]
    assert len(handles) == len(set(handles)) == len(data["handles"])
    assert set(handles) == set(data["handles"]), "missing or substituted handles"
    order = [handles.index(node) for node in data["handles"]]
    raw = np.asarray(native["raw_weights"], dtype=float)[:, order]
    final = np.asarray(native["final_weights"], dtype=float)[:, order]
    assert raw.shape == final.shape == data["fixed"].shape
    assert np.isfinite(raw).all() and np.isfinite(final).all()
    coupled = native["protocol"] == "COUPLED_PARTITION"
    raw_check = oracle.kkt(data, raw, coupled)
    final_check = oracle.kkt(data, final, True)
    assert raw_check["boundary"] <= 1e-9 and raw_check["bounds"] <= 1e-9, raw_check
    for field in ("stationarity", "projected_stationarity", "complementarity"):
        assert raw_check[field] <= 1e-8, (field, raw_check)
    assert final_check["boundary"] <= 1e-9 and final_check["bounds"] <= 1e-9
    assert final_check["partition"] <= 1e-9, final_check
    if coupled:
        assert raw_check["partition"] <= 1e-9
        assert np.array_equal(raw, final), "coupled weights were postprocessed"
        assert native["normalization_change"] == 0
    else:
        assert np.min(raw.sum(axis=1)) > 0
        assert np.max(np.abs(final - raw / raw.sum(axis=1)[:, None])) <= 2e-15
    raw_error = float(np.max(np.abs(raw - np.asarray(expected["raw_weights"]))))
    final_error = float(np.max(np.abs(final - np.asarray(expected["final_weights"]))))
    assert raw_error <= 1e-6 and final_error <= 1e-6, (raw_error, final_error)
    assert abs(native["energy_scale"] / data["scale"] - 1) <= 1e-12
    for named, checked in (("raw_assessment", raw_check), ("final_assessment", final_check)):
        assessment = native[named]
        assert assessment["available"]
        assert abs(assessment["energy"] - checked["energy"]) <= 1e-10 * max(1, checked["energy"])
        for field in ("boundary", "bounds", "partition"):
            assert abs(assessment[field] - checked[field]) <= 1e-12
    attempts = native["attempts"]
    assert len(attempts) == (1 if coupled else len(handles)), "actual call count mismatch"
    assert [x["handle"] for x in attempts] == ([-1] if coupled else list(range(len(handles))))
    assert all(x["candidate_available"] for x in attempts)
    assert all(x["run_invoked"] and x["pass_model_status"] == 0 for x in attempts)
    assert all(isinstance(x["qp_iterations"], int) and x["qp_iterations"] >= 0 for x in attempts)
    if not coupled and native["cells_per_axis"] == 4:
        assert raw_check["partition"] > .09
        assert native["normalization_change"] > .03
        assert final_check["projected_stationarity"] > .01
        assert not native["final_assessment"]["optimality_passed"], "normalized != coupled optimum"
    assert native["success"], "actual run did not validate"
    for field in ("config_sha256", "output_sha256", "fem_output_sha256"):
        assert len(native[field]) == 64 and all(c in "0123456789abcdef" for c in native[field])
    return dict(scope="INDEPENDENT_NATIVE_SYNTHETIC_CHECK_NOT_PAPER_REPRODUCTION",
                cells_per_axis=native["cells_per_axis"], rule=native["rule"], protocol=native["protocol"],
                raw_assessment=raw_check, final_coupled_assessment=final_check,
                max_raw_reference_difference=raw_error, max_final_reference_difference=final_error,
                actual_backend_calls=len(attempts), native_output_sha256=native["output_sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("native")
    parser.add_argument("reference")
    parser.add_argument("--repeat")
    args = parser.parse_args()
    native, reference = load(args.native), load(args.reference)
    result = check(native, reference)
    if args.repeat:
        repeated = load(args.repeat)
        check(repeated, reference)
        assert native == repeated, "scientific output/attempt sequence changed"
        result["deterministic_repeat"] = True
    result["inputs"] = {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                        for p in (args.native, args.reference, *((args.repeat,) if args.repeat else ()))}
    print(json.dumps(result, sort_keys=True, allow_nan=False))
