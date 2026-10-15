from __future__ import annotations

import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "experiments/spin_it_reproduction/strict/compare_faster_repeat.py"
SPEC = importlib.util.spec_from_file_location("faster_repeat_comparison", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise AssertionError("cannot load Faster repeat comparison module")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _result(*, wall_seconds: float, event_seconds: float, objective: float = 1.25):
    return {
        "classification": "RECONSTRUCTED_METHOD_VALIDATION",
        "status": "CONVERGED",
        "state": "LOCKED_PROBLEM_SEALED_AND_INDEPENDENT_AUDITS_PASSED",
        "locked_problem_digest": "a" * 64,
        "function_evaluations": 2,
        "restart_count": 1,
        "trace_hash": "b" * 64,
        "wall_seconds": wall_seconds,
        "events": [{
            "type": "evaluation",
            "evaluation_payload": {
                "objective": objective,
                "evaluator_wall_seconds": event_seconds,
                "residual_evaluation_seconds": event_seconds / 2,
                "derivative_evaluation_seconds": event_seconds / 4,
            },
        }],
        "trace": [{"evaluation": 1, "objective": objective}],
        "evaluation_map_binding_audit": {"status": "CONSISTENT", "bound_evaluation_count": 1},
        "external_terminal_geometry_audit": {
            "status": "PASSED",
            "geometry_verified": True,
            "material_domain_status": "VERIFIED_SUBSET",
            "terminal_material_sha256": "c" * 64,
        },
        "sealed_artifact": {
            "scientific_digest": "d" * 64,
            "problem_digest": "a" * 64,
        },
    }


def _evidence(material: Path, run_record: Path, *, event_objective: float = 1.25,
              source_hash: str = "e" * 64):
    return {
        "begin": {"observation": {
            "source": {"content_sha256": source_hash},
            "runtime": {"native_modules": {
                "_make_it_stand": {"path": "/build/native.so", "sha256": "f" * 64, "size": 10},
            }},
            "dependencies": {"roots": {"cherchi_reference": {"sha256": "a" * 64}}},
        }},
        "scientific": {
            "events": [{"objective": event_objective, "moments": [1.0, 0.0]}],
            "outputs": [
                {"path": str(material), "sha256": MODULE.sha256_file(material), "size": material.stat().st_size},
                {"path": str(run_record), "sha256": MODULE.sha256_file(run_record), "size": run_record.stat().st_size},
            ],
        }
    }


def _write_result(path: Path, record: dict):
    path.write_text(json.dumps(record, sort_keys=True))


class FasterRepeatComparisonTests(unittest.TestCase):
    def test_checkpoint_runtime_context_is_ignored_but_scientific_state_is_bound(self):
        with tempfile.TemporaryDirectory() as temporary:
            roots = [Path(temporary) / name for name in ("first", "repeat")]
            records = [_result(wall_seconds=10.0, event_seconds=.1),
                       _result(wall_seconds=20.0, event_seconds=.9)]
            evidence = []
            for root, record in zip(roots, records):
                root.mkdir()
                # Native hashes include timing; compare their bound event
                # content, not timing-sensitive hashes across processes.
                record["event_hash"] = str(record["wall_seconds"])
                material = root / "material.obj"
                material.write_bytes(b"unchanged scientific material")
                record["evaluations"] = [{"evaluation_index": 0,
                    "oriented_output_path": str(material), "moment_map": {"raw_moments": [1.,0.,0.]}}]
                _write_result(root / "result.json", record)
                observed = _evidence(material, root / "result.json")
                observed["begin"]["observation"]["source"]["git"] = {
                    "diagnostic_elapsed": record["wall_seconds"]}
                observed["begin"]["integrity_sha256"] = hashlib.sha256(json.dumps(
                    observed["begin"],sort_keys=True,separators=(",",":"),
                    ensure_ascii=True,allow_nan=False).encode("ascii")).hexdigest()
                record["evaluation_timings"] = [{"evaluation_index":0,
                    "wall_seconds":record["wall_seconds"],"status":"DEFINED"}]
                record["bridge_scientific_digest"] = hashlib.sha256(
                    json.dumps(record,sort_keys=True).encode()).hexdigest()
                _write_result(root / "result.json", record)
                observed["scientific"]["outputs"][1].update(
                    sha256=MODULE.sha256_file(root/"result.json"),
                    size=(root/"result.json").stat().st_size)
                payload = {
                    "solver_status": "CONVERGED", "trace": record["trace"],
                    "evaluations": record["evaluations"],
                    "evaluation_map_binding_audit": record["evaluation_map_binding_audit"],
                    "design_domain_vertices": [[0.,0.,0.],[1.,0.,0.]],
                    "design_domain_faces": [],
                    "run_finalization": {
                        "schema": "faster-run-finalization-context-v1",
                        "output_directory": str(root), "solver_wall_seconds": record["wall_seconds"],
                        "native_result": {"status":"CONVERGED", "trace":record["trace"],
                            "events":record["events"], "event_hash":record["event_hash"],
                            "evaluation_timings":record["evaluation_timings"],
                            "function_evaluations":2,"restart_count":1},
                        "bridge_scientific_digest":record["bridge_scientific_digest"],
                        "begin":observed["begin"],
                        "classification":record["classification"],
                        "locked_problem_digest":record["locked_problem_digest"],
                        "expected_artifact_digest":record["sealed_artifact"]["scientific_digest"],
                    },
                }
                checkpoint = {"schema":"faster-terminal-audit-checkpoint-v1",
                    "status":"PENDING_TERMINAL_AUDIT", "complete_run_evidence":False,
                    "checker_identity":[str(MODULE_PATH)], "dependencies":[
                        {"path":str(material),"sha256":MODULE.sha256_file(material)},
                        {"path":str(MODULE_PATH),"sha256":MODULE.sha256_file(MODULE_PATH)}],
                    "payload":payload}
                checkpoint["checkpoint_digest"] = hashlib.sha256(json.dumps(
                    checkpoint,sort_keys=True,separators=(",",":"),ensure_ascii=True,
                    allow_nan=False).encode("ascii")).hexdigest()
                path = root / "terminal-audit-checkpoint.json"
                path.write_text(json.dumps(checkpoint,sort_keys=True))
                observed["scientific"]["outputs"].append({"path":str(path),
                    "sha256":MODULE.sha256_file(path),"size":path.stat().st_size})
                evidence.append(observed)
            comparison = MODULE.compare_repeat_records(*records, *evidence,
                first_output_root=roots[0],repeat_output_root=roots[1])
            self.assertEqual(comparison["status"],"SCIENTIFIC_CORE_IDENTICAL")
            self.assertTrue(comparison["checks"]["checkpoint_scientific_states_equal_and_bound"])
            # Invalid checkpoint structures must fail closed, not crash or pass
            # merely because a caller recomputed the file and checkpoint hashes.
            path = roots[1] / "terminal-audit-checkpoint.json"
            original = path.read_text()
            for invalid in ("context_list", "native_list", "empty_checker_identity",
                            "duplicate_dependency", "stale_digest", "changed_dependency",
                            "begin_integrity", "bridge_digest", "timing_identity"):
                with self.subTest(invalid=invalid):
                    checkpoint = json.loads(original)
                    if invalid == "context_list":
                        checkpoint["payload"]["run_finalization"] = []
                    elif invalid == "native_list":
                        checkpoint["payload"]["run_finalization"]["native_result"] = []
                    elif invalid == "empty_checker_identity":
                        checkpoint["checker_identity"] = []
                    elif invalid == "duplicate_dependency":
                        checkpoint["dependencies"].append(dict(checkpoint["dependencies"][0]))
                    elif invalid == "changed_dependency":
                        checkpoint["dependencies"][0]["sha256"] = "0" * 64
                    elif invalid == "begin_integrity":
                        checkpoint["payload"]["run_finalization"]["begin"]["integrity_sha256"] = "0" * 64
                    elif invalid == "bridge_digest":
                        checkpoint["payload"]["run_finalization"]["bridge_scientific_digest"] = "0" * 64
                    elif invalid == "timing_identity":
                        checkpoint["payload"]["run_finalization"]["native_result"]["evaluation_timings"][0]["evaluation_index"] = 9
                    else:
                        checkpoint["payload"]["design_domain_vertices"][0][0] = 3.
                    if invalid != "stale_digest":
                        checkpoint.pop("checkpoint_digest")
                        checkpoint["checkpoint_digest"] = hashlib.sha256(json.dumps(
                            checkpoint,sort_keys=True,separators=(",",":"),ensure_ascii=True,
                            allow_nan=False).encode("ascii")).hexdigest()
                    path.write_text(json.dumps(checkpoint,sort_keys=True))
                    evidence[1]["scientific"]["outputs"][-1].update(
                        sha256=MODULE.sha256_file(path),size=path.stat().st_size)
                    if invalid in ("empty_checker_identity", "duplicate_dependency"):
                        with self.assertRaises(ValueError):
                            MODULE._checkpoint_scientific_state(roots[1],records[1],evidence[1])
                    rejected = MODULE.compare_repeat_records(*records,*evidence,
                        first_output_root=roots[0],repeat_output_root=roots[1])
                    self.assertEqual(rejected["status"],"MISMATCH")
                    self.assertFalse(rejected["checks"]["checkpoint_scientific_states_equal_and_bound"])
            path.write_text(original)
            evidence[1]["scientific"]["outputs"][-1].update(
                sha256=MODULE.sha256_file(path),size=path.stat().st_size)
            # Rehashing a changed checkpoint must not hide a design-domain change.
            path = roots[1] / "terminal-audit-checkpoint.json"
            checkpoint = json.loads(path.read_text())
            checkpoint["payload"]["design_domain_vertices"][1][0] = 2.
            checkpoint.pop("checkpoint_digest")
            checkpoint["checkpoint_digest"] = hashlib.sha256(json.dumps(
                checkpoint,sort_keys=True,separators=(",",":"),ensure_ascii=True,
                allow_nan=False).encode("ascii")).hexdigest()
            path.write_text(json.dumps(checkpoint,sort_keys=True))
            evidence[1]["scientific"]["outputs"][-1].update(
                sha256=MODULE.sha256_file(path),size=path.stat().st_size)
            changed = MODULE.compare_repeat_records(*records,*evidence,
                first_output_root=roots[0],repeat_output_root=roots[1])
            self.assertEqual(changed["status"],"MISMATCH")
            self.assertFalse(changed["checks"]["checkpoint_scientific_states_equal_and_bound"])

    def test_volatile_elapsed_times_do_not_break_identical_scientific_run_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_root = root / "first" / "solver-output"
            repeat_root = root / "repeat" / "solver-output"
            first_root.mkdir(parents=True)
            repeat_root.mkdir(parents=True)
            first_material = first_root / "material.obj"
            repeat_material = repeat_root / "material.obj"
            first_result_file = first_root / "result.json"
            repeat_result_file = repeat_root / "result.json"
            first_material.write_bytes(b"same sealed scientific mesh")
            repeat_material.write_bytes(b"same sealed scientific mesh")
            first_result = _result(wall_seconds=10.0, event_seconds=0.1)
            repeat_result = _result(wall_seconds=20.0, event_seconds=0.9)
            _write_result(first_result_file, first_result)
            _write_result(repeat_result_file, repeat_result)

            comparison = MODULE.compare_repeat_records(
                first_result,
                repeat_result,
                _evidence(first_material, first_result_file),
                _evidence(repeat_material, repeat_result_file),
                first_output_root=first_root,
                repeat_output_root=repeat_root,
            )

        self.assertEqual(comparison["status"], "SCIENTIFIC_CORE_IDENTICAL")
        self.assertTrue(comparison["checks"]["solver_events_equal_without_timing"])
        self.assertTrue(comparison["checks"]["sealed_geometry_and_moments_identical"])
        self.assertEqual(comparison["outputs"]["identical_scientific_file_count"], 1)
        self.assertEqual(comparison["outputs"]["different_runtime_metadata_files"], ["result.json"])

    def test_scientific_event_change_fails_repeat_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_root = root / "first" / "solver-output"
            repeat_root = root / "repeat" / "solver-output"
            first_root.mkdir(parents=True)
            repeat_root.mkdir(parents=True)
            first_material = first_root / "material.obj"
            repeat_material = repeat_root / "material.obj"
            first_result_file = first_root / "result.json"
            repeat_result_file = repeat_root / "result.json"
            first_material.write_bytes(b"same sealed scientific mesh")
            repeat_material.write_bytes(b"same sealed scientific mesh")
            first_result = _result(wall_seconds=10.0, event_seconds=0.1)
            repeat_result = _result(wall_seconds=20.0, event_seconds=0.9, objective=1.5)
            _write_result(first_result_file, first_result)
            _write_result(repeat_result_file, repeat_result)

            comparison = MODULE.compare_repeat_records(
                first_result,
                repeat_result,
                _evidence(first_material, first_result_file),
                _evidence(repeat_material, repeat_result_file, event_objective=1.5),
                first_output_root=first_root,
                repeat_output_root=repeat_root,
            )

        self.assertEqual(comparison["status"], "MISMATCH")
        self.assertFalse(comparison["checks"]["solver_events_equal_without_timing"])

    def test_scientific_geometry_byte_change_fails_repeat_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_root = root / "first" / "solver-output"
            repeat_root = root / "repeat" / "solver-output"
            first_root.mkdir(parents=True)
            repeat_root.mkdir(parents=True)
            first_material = first_root / "material.obj"
            repeat_material = repeat_root / "material.obj"
            first_result_file = first_root / "result.json"
            repeat_result_file = repeat_root / "result.json"
            first_material.write_bytes(b"first mesh")
            repeat_material.write_bytes(b"different mesh")
            first_result = _result(wall_seconds=10.0, event_seconds=0.1)
            repeat_result = _result(wall_seconds=20.0, event_seconds=0.9)
            _write_result(first_result_file, first_result)
            _write_result(repeat_result_file, repeat_result)

            comparison = MODULE.compare_repeat_records(
                first_result,
                repeat_result,
                _evidence(first_material, first_result_file),
                _evidence(repeat_material, repeat_result_file),
                first_output_root=first_root,
                repeat_output_root=repeat_root,
            )

        self.assertEqual(comparison["status"], "MISMATCH")
        self.assertFalse(comparison["checks"]["scientific_output_bytes_equal"])

    def test_mismatched_effective_source_rejects_repeat_claim(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_root = root / "first" / "solver-output"
            repeat_root = root / "repeat" / "solver-output"
            first_root.mkdir(parents=True)
            repeat_root.mkdir(parents=True)
            for output_root in (first_root, repeat_root):
                (output_root / "material.obj").write_bytes(b"same scientific material")
            first_result = _result(wall_seconds=10.0, event_seconds=0.1)
            repeat_result = _result(wall_seconds=20.0, event_seconds=0.9)
            _write_result(first_root / "result.json", first_result)
            _write_result(repeat_root / "result.json", repeat_result)
            comparison = MODULE.compare_repeat_records(
                first_result,
                repeat_result,
                _evidence(first_root / "material.obj", first_root / "result.json"),
                _evidence(repeat_root / "material.obj", repeat_root / "result.json",
                          source_hash="d" * 64),
                first_output_root=first_root,
                repeat_output_root=repeat_root,
            )

        self.assertEqual(comparison["status"], "MISMATCH")
        self.assertFalse(comparison["checks"]["source_binary_and_dependencies_match"])

    def test_comparison_digest_binds_the_verified_output_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            comparisons = []
            for case, material_bytes in (("alpha", b"first valid material"),
                                         ("beta", b"second valid material")):
                output_roots = [root / case / label / "solver-output"
                                for label in ("first", "repeat")]
                first_result = _result(wall_seconds=10.0, event_seconds=0.1)
                repeat_result = _result(wall_seconds=20.0, event_seconds=0.9)
                for output_root in output_roots:
                    output_root.mkdir(parents=True)
                    (output_root / "material.obj").write_bytes(material_bytes)
                _write_result(output_roots[0] / "result.json", first_result)
                _write_result(output_roots[1] / "result.json", repeat_result)
                comparison = MODULE.compare_repeat_records(
                    first_result,
                    repeat_result,
                    _evidence(output_roots[0] / "material.obj", output_roots[0] / "result.json"),
                    _evidence(output_roots[1] / "material.obj", output_roots[1] / "result.json"),
                    first_output_root=output_roots[0],
                    repeat_output_root=output_roots[1],
                )
                self.assertEqual(comparison["status"], "SCIENTIFIC_CORE_IDENTICAL")
                comparisons.append(comparison)

        self.assertNotEqual(comparisons[0]["comparison_digest"],
                            comparisons[1]["comparison_digest"])

    def test_result_object_must_match_evidence_bound_result_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_root = root / "first" / "solver-output"
            repeat_root = root / "repeat" / "solver-output"
            first_root.mkdir(parents=True)
            repeat_root.mkdir(parents=True)
            first_result = _result(wall_seconds=10.0, event_seconds=0.1)
            repeat_result = _result(wall_seconds=20.0, event_seconds=0.9)
            first_file = first_root / "result.json"
            repeat_file = repeat_root / "result.json"
            first_file.write_text(json.dumps(first_result))
            altered_on_disk = {**repeat_result, "status": "FAILED"}
            repeat_file.write_text(json.dumps(altered_on_disk))
            first_material = first_root / "material.obj"
            repeat_material = repeat_root / "material.obj"
            first_material.write_bytes(b"same material")
            repeat_material.write_bytes(b"same material")
            comparison = MODULE.compare_repeat_records(
                first_result, repeat_result,
                _evidence(first_material, first_file),
                _evidence(repeat_material, repeat_file),
                first_output_root=first_root,
                repeat_output_root=repeat_root,
            )

        self.assertEqual(comparison["status"], "MISMATCH")
        self.assertFalse(comparison["checks"]["result_records_match_verified_files"])


import copy

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),
        ensure_ascii=True,allow_nan=False).encode('ascii')).hexdigest()

def fixture(root,seconds):
    root.mkdir();material=root/'material.obj';material.write_bytes(b'same-scientific-material')
    result=_result(wall_seconds=seconds,event_seconds=seconds/100)
    result.update(evaluations=[],event_hash=str(seconds),
                  evaluation_timings=[{'evaluation_index':0,'wall_seconds':seconds,'status':'DEFINED'}])
    _write_result(root/'result.json',result)
    evidence=_evidence(material,root/'result.json')
    native={key:copy.deepcopy(result[key]) for key in
        ('status','trace','events','event_hash','evaluation_timings','function_evaluations','restart_count')}
    bridge={'operation':'SOLVE_CHERCHI_EXTERNAL_GEOMETRY','native_run':copy.deepcopy(native),
            'status':'CONVERGED','seal_error':None,'problem_digest':result['locked_problem_digest'],
            'artifact_scientific_digest':result['sealed_artifact']['scientific_digest']}
    context={'schema':'faster-run-finalization-context-v1','output_directory':str(root),
        'solver_wall_seconds':seconds,'native_result':native,'bridge_identity':bridge,
        'bridge_status':'CONVERGED','bridge_seal_error':None,'begin':evidence['begin'],
        'classification':result['classification'],'locked_problem_digest':result['locked_problem_digest'],
        'expected_artifact_digest':result['sealed_artifact']['scientific_digest']}
    checkpoint={'schema':'faster-terminal-audit-checkpoint-v1','status':'PENDING_TERMINAL_AUDIT',
        'complete_run_evidence':False,'checker_identity':[str(MODULE_PATH)],
        'dependencies':[{'path':str(MODULE_PATH),'sha256':MODULE.sha256_file(MODULE_PATH)}],
        'payload':{'solver_status':'CONVERGED','trace':result['trace'],'evaluations':[],
                   'evaluation_map_binding_audit':result['evaluation_map_binding_audit'],
                   'run_finalization':context}}
    return result,evidence,checkpoint

def persist(root,result,evidence,checkpoint):
    context=checkpoint['payload']['run_finalization']
    bridge=context['bridge_identity']
    context['bridge_scientific_digest']=result['bridge_scientific_digest']=digest(bridge)
    result['bridge_identity_attachment']={'schema':'cherchi-bridge-identity-attachment-v1',
        'status':'BOUND_IDENTITY','identity':copy.deepcopy(bridge),'digest':digest(bridge)}
    checkpoint.pop('checkpoint_digest',None);checkpoint['checkpoint_digest']=digest(checkpoint)
    _write_result(root/'result.json',result)
    _write_result(root/'terminal-audit-checkpoint.json',checkpoint)
    evidence['scientific']['outputs']=[{'path':str(root/name),'sha256':MODULE.sha256_file(root/name),
        'size':(root/name).stat().st_size} for name in ('material.obj','result.json','terminal-audit-checkpoint.json')]

class BridgeCheckpointTests(unittest.TestCase):
    def test_duplicate_native_timing_copy_is_normalized_only_after_binding(self):
        with tempfile.TemporaryDirectory() as temp:
            roots=[Path(temp)/name for name in ('first','repeat')]
            rows=[fixture(root,t) for root,t in zip(roots,(10.,20.))]
            for root,row in zip(roots,rows):persist(root,*row)
            result=MODULE.compare_repeat_records(rows[0][0],rows[1][0],rows[0][1],rows[1][1],
                first_output_root=roots[0],repeat_output_root=roots[1])
            self.assertEqual(result['status'],'SCIENTIFIC_CORE_IDENTICAL')

    def test_rehashed_inconsistent_bridge_copy_is_not_accepted(self):
        for mutation in ('native_count','timing_index','problem','missing_native'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as temp:
                root=Path(temp)/'run';result,evidence,checkpoint=fixture(root,10.)
                bridge=checkpoint['payload']['run_finalization']['bridge_identity']
                if mutation=='native_count':bridge['native_run']['function_evaluations']=99
                elif mutation=='timing_index':bridge['native_run']['evaluation_timings'][0]['evaluation_index']=99
                elif mutation=='problem':bridge['problem_digest']='0'*64
                else:bridge.pop('native_run')
                persist(root,result,evidence,checkpoint)
                with self.assertRaises(ValueError):MODULE._checkpoint_scientific_state(root,result,evidence)

    def test_unknown_scientific_bridge_field_is_not_discarded(self):
        with tempfile.TemporaryDirectory() as temp:
            roots=[Path(temp)/name for name in ('first','repeat')]
            rows=[fixture(root,10.) for root in roots]
            for index,(root,row) in enumerate(zip(roots,rows)):
                row[2]['payload']['run_finalization']['bridge_identity']['future_scientific_field']=index
                persist(root,*row)
            result=MODULE.compare_repeat_records(rows[0][0],rows[1][0],rows[0][1],rows[1][1],
                first_output_root=roots[0],repeat_output_root=roots[1])
            self.assertEqual(result['status'],'MISMATCH')
            self.assertFalse(result['checks']['checkpoint_scientific_states_equal_and_bound'])

class BridgePresenceTests(unittest.TestCase):
    def test_one_sided_identity_removal_cannot_be_downgraded_to_legacy(self):
        for missing in ("context", "result"):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as temp:
                root=Path(temp)/"run"
                result,evidence,checkpoint=fixture(root,10.)
                persist(root,result,evidence,checkpoint)
                if missing=="context":
                    checkpoint["payload"]["run_finalization"].pop("bridge_identity")
                else:
                    result.pop("bridge_identity_attachment")
                checkpoint.pop("checkpoint_digest")
                checkpoint["checkpoint_digest"]=digest(checkpoint)
                _write_result(root/"result.json",result)
                _write_result(root/"terminal-audit-checkpoint.json",checkpoint)
                evidence["scientific"]["outputs"]=[
                    {"path":str(root/name),"sha256":MODULE.sha256_file(root/name),
                     "size":(root/name).stat().st_size}
                    for name in ("material.obj","result.json","terminal-audit-checkpoint.json")]
                with self.assertRaises(ValueError):
                    MODULE._checkpoint_scientific_state(root,result,evidence)

if __name__ == "__main__":
    unittest.main()
