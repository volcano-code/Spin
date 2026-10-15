from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import unittest
import tempfile
import hashlib
import subprocess
import sys
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "experiments/spin_it_reproduction/faster/run_paper_figures.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("spin_it_faster_figures", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Faster paper-figure runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpinItFasterPaperFigureTests(unittest.TestCase):
    def bound_journal_fixture(self, directory):
        root = Path(directory)
        output = root/'solver-output'
        output.mkdir()
        journal = output/'solver-events.jsonl'
        journal.write_text(''.join(json.dumps(row)+'\n' for row in self.journal_rows()))
        manifest = root/'problem.json'
        manifest.write_text(json.dumps({'schema': 'faster-cherchi-problem-manifest-v1',
                                        'target': {'mode': 'YO_YO'}}))
        config = root/'config.json'
        config.write_text(json.dumps({
            'argv': [sys.executable, 'runner.py', '--problem-manifest', str(manifest),
                     '--output', str(output)],
            'identity': {str(manifest): {'path': str(manifest),
                'sha256': hashlib.sha256(manifest.read_bytes()).hexdigest(),
                'size': manifest.stat().st_size}},
        }))
        return journal, manifest, config

    def test_figure10_run_config_derives_yoyo_group_without_caller_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            journal, manifest, config = self.bound_journal_fixture(directory)
            payload = load_runner().build_payload_from_solver_journal(
                journal, case_id='method-only', run_config=config)
            self.assertEqual(payload['figure10']['groups']['yo_yos'], [(0.1, 1e-13)])
            self.assertEqual(payload['figure10']['groups']['spinning_tops'], [])
            binding = payload['figure10']['solver_journal']['application_binding']
            self.assertEqual(binding['mode'], 'YO_YO')
            self.assertEqual(binding['manifest_sha256'], hashlib.sha256(manifest.read_bytes()).hexdigest())
            self.assertEqual(binding['run_config_sha256'], hashlib.sha256(config.read_bytes()).hexdigest())
            self.assertFalse(payload['author_exact'])
            self.assertEqual(payload['figure10']['solver_journal']['input_identity_attestation'],
                             'NOT_VERIFIED_BY_FIGURE_CHECKER')

    def test_figure10_run_config_rejects_caller_mislabeling_yoyo_as_top(self):
        with tempfile.TemporaryDirectory() as directory:
            journal, _, config = self.bound_journal_fixture(directory)
            with self.assertRaisesRegex(ValueError, 'mode'):
                load_runner().build_payload_from_solver_journal(
                    journal, case_id='method-only', mode='SPINNING_TOP', run_config=config)

    def test_figure10_run_config_rejects_changed_manifest_and_foreign_journal(self):
        with tempfile.TemporaryDirectory() as directory:
            journal, manifest, config = self.bound_journal_fixture(directory)
            foreign = Path(directory)/'foreign.jsonl'
            foreign.write_bytes(journal.read_bytes())
            with self.assertRaisesRegex(ValueError, 'journal'):
                load_runner().build_payload_from_solver_journal(
                    foreign, case_id='method-only', run_config=config)
            manifest.write_text(json.dumps({'schema': 'faster-cherchi-problem-manifest-v1',
                                            'target': {'mode': 'SPINNING_TOP'}}))
            with self.assertRaisesRegex(ValueError, 'identity'):
                load_runner().build_payload_from_solver_journal(
                    journal, case_id='method-only', run_config=config)

    def test_figure10_run_config_rejects_missing_or_ambiguous_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            journal, manifest, config = self.bound_journal_fixture(directory)
            original = json.loads(config.read_text())
            invalid = []
            for option in ('--problem-manifest', '--output'):
                removed = json.loads(json.dumps(original))
                index = removed['argv'].index(option)
                del removed['argv'][index:index+2]
                invalid.append(removed)
                duplicate = json.loads(json.dumps(original))
                duplicate['argv'] += [option, original['argv'][index+1]]
                invalid.append(duplicate)
            unpinned = json.loads(json.dumps(original))
            unpinned['identity'] = {}
            invalid.append(unpinned)
            for bad in invalid:
                config.write_text(json.dumps(bad))
                with self.subTest(bad=bad), self.assertRaises(ValueError):
                    load_runner().build_payload_from_solver_journal(
                        journal, case_id='method-only', run_config=config)

    def test_figure10_run_config_cli_uses_actual_mode_and_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            journal, _, config = self.bound_journal_fixture(directory)
            output = Path(directory)/'figure.json'
            command = [sys.executable, str(RUNNER), '--solver-journal', str(journal),
                       '--journal-case-id', 'method-only', '--run-config', str(config),
                       '--output', str(output)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(output.read_text())
            self.assertEqual(payload['figure10']['groups']['yo_yos'], [[0.1, 1e-13]])
            mismatch = subprocess.run(command+['--journal-mode', 'SPINNING_TOP'],
                                      capture_output=True, text=True, timeout=30)
            self.assertNotEqual(mismatch.returncode, 0)
            self.assertIn('mode', mismatch.stderr)

    def test_figure10_run_config_rejects_argparse_abbreviation_overrides(self):
        with tempfile.TemporaryDirectory() as directory:
            journal, _, config = self.bound_journal_fixture(directory)
            original = json.loads(config.read_text())
            # The actual solver's argparse accepts these abbreviated options
            # and the later value wins. Never bind a different literal option.
            for extra in (['--problem', '/different-top.json'],
                          ['--out', '/different-output'],
                          ['--problem=/different-top.json'],
                          ['--out=/different-output']):
                config.write_text(json.dumps({**original, 'argv': original['argv']+extra}))
                with self.subTest(extra=extra), self.assertRaises(ValueError):
                    load_runner().build_payload_from_solver_journal(
                        journal, case_id='method-only', run_config=config)

    def test_figure10_run_config_cannot_be_silently_ignored_without_journal(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, config = self.bound_journal_fixture(directory)
            output = Path(directory)/'figure.json'
            result = subprocess.run([sys.executable, str(RUNNER), '--run-config', str(config),
                                     '--output', str(output)],
                                    capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())

    def journal_rows(self):
        def evaluation(index, outcome, raw=None):
            return {'record_type': 'solver_event', 'record': {
                'type': 'evaluation', 'evaluation_index': index, 'restart': 0,
                'iteration': 0 if index == 1 else 1, 'outcome': outcome,
                'reason': 'evaluation_exception' if outcome == 'failed' else 'residual_decreased',
                'raw_equation_residual': raw, 'equation_profile': 'DUAL_MONOMIAL_RAW_V1',
                'multiplier_mu': 2., 'physical_relative_residual': 1e-9}}
        return [
            {'schema': 'faster-streamed-solver-events-v1', 'inputs_sha256': 'a'*64},
            evaluation(1, 'accepted', 0.1), evaluation(2, 'failed'),
            evaluation(3, 'accepted', 1e-13),
            {'record_type': 'solver_event', 'record': {
                'type': 'attempt_terminal', 'restart': 0, 'evaluation_index': 3,
                'outcome': 'terminated', 'reason': 'converged'}},
        ]

    def journal_payload(self, rows, *, trailing_newline=True):
        runner = load_runner()
        checker = getattr(runner, 'build_payload_from_solver_journal', None)
        self.assertTrue(callable(checker), 'actual solver journal Figure adapter is missing')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'solver-events.jsonl'
            content = '\n'.join(json.dumps(row) for row in rows)
            if trailing_newline:
                content += '\n'
            path.write_text(content)
            return checker(path, case_id='method-not-author-asset', mode='YO_YO'), content

    def test_figure10_actual_journal_retains_failed_evaluation_and_source_identity(self):
        payload, content = self.journal_payload(self.journal_rows())
        source = payload['figure10']['solver_journal']
        self.assertEqual(source['sha256'], hashlib.sha256(content.encode()).hexdigest())
        self.assertEqual(source['inputs_sha256'], 'a'*64)
        self.assertEqual(source['evaluations'], 3)
        self.assertEqual(source['failed_evaluations'], 1)
        self.assertEqual(source['terminal_attempts'], 1)
        self.assertEqual(source['attempts'][0]['numerical_acceptance']['passed'], True)
        self.assertEqual(payload['figure10']['groups']['yo_yos'], [(0.1, 1e-13)])
        self.assertFalse(payload['author_exact'])
        self.assertEqual(source['scope'], 'JOURNAL_NUMERICS_NOT_MATERIAL_OR_ORIGINAL_PAPER_ACCEPTANCE')

    def test_figure10_accepts_typed_diagnostic_journal_without_counting_it_as_convergence(self):
        rows = self.journal_rows()
        rows[0]['schema'] = 'faster-streamed-solver-events-v2'
        rows[2]['record']['newton_system_payload'] = {
            'linear_backward_error': {'$nonfinite_float': '+inf'}}
        payload, _ = self.journal_payload(rows)
        source = payload['figure10']['solver_journal']
        self.assertEqual(source['failed_evaluations'], 1)
        self.assertEqual(source['evaluations'], 3)
        self.assertEqual(payload['figure10']['groups']['yo_yos'], [(0.1, 1e-13)])
        self.assertFalse(payload['author_exact'])

    def test_figure10_unterminated_attempt_cannot_supply_successful_pairs(self):
        payload, _ = self.journal_payload(self.journal_rows()[:-1])
        source = payload['figure10']['solver_journal']
        self.assertEqual(source['attempts'][0]['status'], 'INCOMPLETE')
        self.assertEqual(source['terminal_attempts'], 0)
        self.assertEqual(payload['figure10']['groups']['yo_yos'], [])

    def test_figure10_journal_rejects_dropped_evaluation_and_truncated_line(self):
        rows = self.journal_rows()
        for invalid in (rows[:2]+rows[3:], rows[:2]+rows[1:]):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.journal_payload(invalid)
        with self.assertRaises(ValueError):
            self.journal_payload(rows, trailing_newline=False)

    def test_figure10_journal_does_not_approve_negative_multiplier_terminal(self):
        rows = self.journal_rows()
        rows[3]['record']['multiplier_mu'] = -2.
        payload, _ = self.journal_payload(rows)
        check = payload['figure10']['solver_journal']['attempts'][0]['numerical_acceptance']
        self.assertFalse(check['passed'])
        self.assertIn('POSITIVE_MULTIPLIER_REQUIRED', check['reasons'])
        self.assertEqual(payload['figure10']['groups']['yo_yos'], [])

    def test_figure10_journal_rejects_missing_restart_attempt(self):
        rows = self.journal_rows()
        rows.append({'record_type': 'solver_event', 'record': {
            'type': 'attempt_terminal', 'restart': 2, 'evaluation_index': 3,
            'outcome': 'terminated', 'reason': 'initialization_failed'}})
        with self.assertRaises(ValueError):
            self.journal_payload(rows)

    def test_figure10_pairs_successful_monotone_residuals_by_application(self) -> None:
        runner = load_runner()
        table = {
            "cases": [
                {
                    "id": "static",
                    "mode": "STATIC_BALANCE",
                    "attempts": [{
                        "status": "CONVERGED",
                        "trace": [
                            {"raw_equation_residual": 1e-1, "normalized_equation_residual": 0.5},
                            {"raw_equation_residual": 1e-2, "normalized_equation_residual": 0.25},
                            {"raw_equation_residual": 1e-4, "normalized_equation_residual": 0.125},
                        ],
                    }],
                },
                {
                    "id": "failed",
                    "mode": "SPINNING_TOP",
                    "attempts": [{"status": "RESOURCE_LIMIT_REACHED", "trace": []}],
                },
            ]
        }
        grouped = runner.quadratic_convergence_pairs(table)
        self.assertEqual(grouped["static_stability"], [(1e-1, 1e-2), (1e-2, 1e-4)])
        self.assertEqual(grouped["spinning_tops"], [])

    def test_figure10_cannot_replace_missing_raw_residual_with_normalized(self) -> None:
        runner = load_runner()
        for bad in (None, True, float("nan"), float("inf"), -1.0):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                runner.quadratic_convergence_pairs({"cases": [{
                    "mode": "YO_YO", "attempts": [{"status": "CONVERGED", "trace": [
                        {"normalized_equation_residual": 1e-5, "raw_equation_residual": bad}
                    ]}]
                }]})

    def test_figure10_order_is_estimated_only_above_noise(self) -> None:
        runner = load_runner()
        table = {"cases": [{"mode": "SPINNING_TOP", "attempts": [{
            "status": "CONVERGED", "trace": [
                {"raw_equation_residual": value}
                for value in (0.1, 0.01, 0.0001, 1e-8, 1e-16)
            ]
        }]}]}
        report = runner.convergence_order_evidence(table, noise_floor=1e-12)
        evidence = report["spinning_tops"][0]
        self.assertEqual(evidence["status"], "ESTIMATED")
        self.assertEqual(evidence["usable_triples"], 2)
        self.assertEqual(evidence["excluded_triples"], 1)
        self.assertTrue(all(abs(value - 2.0) < 1e-12 for value in evidence["orders"]))
        table["cases"][0]["attempts"][0]["trace"] = [
            {"raw_equation_residual": value} for value in (1e-8, 1e-16)
        ]
        insufficient = runner.convergence_order_evidence(table)["spinning_tops"][0]
        self.assertEqual(insufficient["status"], "INSUFFICIENT_EVIDENCE")
        self.assertIsNone(insufficient["median_order"])

    def test_figure10_does_not_join_different_restart_trajectories(self) -> None:
        runner = load_runner()
        table = {"cases": [{"mode": "YO_YO", "attempts": [{
            "status": "CONVERGED", "trace": [
                {"restart": 0, "raw_equation_residual": 0.1},
                {"restart": 0, "raw_equation_residual": 0.01},
                {"restart": 1, "raw_equation_residual": 0.2},
                {"restart": 1, "raw_equation_residual": 0.04},
            ]
        }]}]}
        self.assertEqual(runner.quadratic_convergence_pairs(table)["yo_yos"],
                         [(0.1, 0.01), (0.2, 0.04)])
        order = runner.convergence_order_evidence(table)["yo_yos"][0]
        self.assertEqual(order["status"], "INSUFFICIENT_EVIDENCE")
        self.assertEqual(order["usable_triples"], 0)

    def test_figure10_late_order_does_not_borrow_quadratic_failed_restart(self) -> None:
        runner = load_runner()
        table = {"cases": [{"mode": "SPINNING_TOP", "attempts": [{
            "status": "CONVERGED", "trace": [
                {"restart": restart, "raw_equation_residual": value}
                for restart, values in (
                    (0, (0.1, 0.01, 0.0001, 1e-8)),
                    (1, (0.1, 0.05, 0.025, 0.0125, 0.00625, 1e-13)),
                ) for value in values
            ]
        }]}]}
        evidence = runner.convergence_order_evidence(table)["spinning_tops"][0]
        tail = evidence["terminal_restart_order_evidence"]
        self.assertEqual(tail["restart"], 1)
        self.assertEqual(tail["status"], "ESTIMATED")
        self.assertEqual(tail["usable_triples"], 3)
        self.assertAlmostEqual(tail["median_order"], 1.0)
        self.assertTrue(all(abs(order - 1.0) < 1e-12 for order in tail["orders"]))

    def test_figure10_late_order_retains_only_predeclared_terminal_window(self) -> None:
        runner = load_runner()
        values = (0.1, 0.01, 0.0001, 1e-5, 5e-6, 2.5e-6, 1.25e-6,
                  6.25e-7, 3.125e-7, 1.5625e-7, 1e-13)
        table = {"cases": [{"mode": "YO_YO", "attempts": [{
            "status": "CONVERGED", "trace": [
                {"raw_equation_residual": value} for value in values
            ]
        }]}]}
        tail = runner.convergence_order_evidence(table)["yo_yos"][0][
            "terminal_restart_order_evidence"]
        self.assertEqual(tail["maximum_triples"], 5)
        self.assertEqual(tail["usable_triples"], 5)
        self.assertTrue(all(abs(order - 1.0) < 1e-12 for order in tail["orders"]))

    def test_figure10_rejects_restart_identity_reuse_and_regression(self) -> None:
        runner = load_runner()
        for restarts in ((0, 1, 0), (2, 1)):
            with self.subTest(restarts=restarts), self.assertRaises(ValueError):
                runner.quadratic_convergence_pairs({"cases": [{"mode": "YO_YO",
                    "attempts": [{"status": "CONVERGED", "trace": [
                        {"restart": restart, "raw_equation_residual": 0.1}
                        for restart in restarts]}]}]})

    def test_figure10_late_order_does_not_bridge_a_plateau(self) -> None:
        runner = load_runner()
        table = {"cases": [{"mode": "YO_YO", "attempts": [{
            "status": "CONVERGED", "trace": [
                {"raw_equation_residual": value}
                for value in (0.1, 0.05, 0.025, 0.025, 0.0125, 0.00625)
            ]
        }]}]}
        tail = runner.convergence_order_evidence(table)["yo_yos"][0][
            "terminal_restart_order_evidence"]
        self.assertEqual(tail["status"], "INSUFFICIENT_EVIDENCE")
        self.assertEqual(tail["usable_triples"], 1)

    def test_figure10_fixed_branch_order_does_not_bridge_events_or_missing_identity(self):
        runner = load_runner()
        values = (0.1, 0.01, 0.0001, 1e-8, 1e-16)
        def check(branches):
            table = {"cases": [{"mode": "YO_YO", "attempts": [{
                "status": "CONVERGED", "trace": [
                    {"raw_equation_residual": value, "material_topology_hash": branch}
                    for value, branch in zip(values, branches)
                ]}]}]}
            return runner.convergence_order_evidence(table)["yo_yos"][0]
        same = check(['a'*64]*5)
        self.assertIn('fixed_branch_terminal_order_evidence', same)
        fixed = same['fixed_branch_terminal_order_evidence']
        self.assertEqual(fixed['usable_triples'], 2)
        self.assertEqual(fixed['status'], 'ESTIMATED')
        self.assertAlmostEqual(fixed['median_order'], 2.)
        # The old whole-trajectory diagnostic stays visible, but cannot lend
        # its earlier samples to the branch at the actual endpoint.
        for branches, count in ((['a'*64]+['b'*64]*4, 1),
                                (['a'*64]*4+['b'*64], 0),
                                (['a'*64, None]+['a'*64]*3, 0),
                                ([None]*5, 0)):
            with self.subTest(branches=branches):
                evidence = check(branches)
                self.assertEqual(evidence['status'], 'ESTIMATED')
                fixed = evidence['fixed_branch_terminal_order_evidence']
                self.assertEqual(fixed['status'], 'INSUFFICIENT_EVIDENCE')
                self.assertEqual(fixed['usable_triples'], count)
        for invalid in ('', 'a'*63, 'G'*64, True, {}):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                check(['a'*64]*4+[invalid])

    def test_figure10_journal_retains_actual_material_branch_identity(self):
        rows = self.journal_rows()
        for row in rows:
            record = row.get('record', {})
            if record.get('outcome') == 'accepted':
                record['material_topology_hash'] = 'a'*64
        payload, _ = self.journal_payload(rows)
        trace = payload['figure10']['solver_journal']['attempts'][0]['trace']
        self.assertEqual([row.get('material_topology_hash') for row in trace], ['a'*64]*2)

    def test_pareto_weight_grid_and_front_invariants(self) -> None:
        runner = load_runner()
        weights = runner.pareto_weights(11)
        self.assertEqual(weights[0], 0.0)
        self.assertEqual(weights[-1], 1.0)
        self.assertTrue(all(a < b for a, b in zip(weights, weights[1:])))
        front = runner.validate_pareto_front([
            {"weight": 0.0, "first_objective": 4.0, "center_height": 0.2},
            {"weight": 0.5, "first_objective": 3.0, "center_height": 0.3},
            {"weight": 1.0, "first_objective": 2.0, "center_height": 0.4},
        ])
        self.assertTrue(front["finite"])
        self.assertTrue(front["nondominated"])

    def test_figure13_binary_search_and_distances_match_paper_definition(self) -> None:
        runner = load_runner()
        critical = runner.binary_search_critical_distance(
            lambda distance: distance <= 1.0,
            lower=0.0,
            upper=2.0,
            tolerance=1e-9,
        )
        self.assertLessEqual(abs(critical - 1.0), 1e-9)
        distances = runner.feasibility_distances(critical)
        self.assertEqual(len(distances), 9)
        for exponent, distance in enumerate(distances, 1):
            self.assertAlmostEqual(critical - distance, 10.0 ** (-exponent))

    def test_degenerate_one_dimensional_examples_are_exact(self) -> None:
        runner = load_runner()
        examples = runner.degenerate_one_dimensional_examples()
        self.assertEqual(examples[0]["optimizer_cardinality"], "UNCOUNTABLY_INFINITE")
        self.assertEqual(examples[1]["optimizer_cardinality"], "ONE")
        self.assertEqual(examples[1]["unique_fill"], 1.0)
        self.assertTrue(all(math.isclose(row["gradient_at_optimum"], 0.0) for row in examples))

    def test_figure13_missing_continuous_oracle_does_not_run_newton(self) -> None:
        runner = load_runner()
        with patch("make_it_stand.spin.solve_spin_it_faster", side_effect=AssertionError("no proof")):
            result = runner.run_tilted_cube_feasibility(tolerance=0.05)
        self.assertEqual(result["feasibility_status"], "UNKNOWN")
        self.assertIsNone(result["critical_distance"])
        self.assertEqual(result["rows"], [])
        self.assertIn("CONTINUOUS_FEASIBILITY_ORACLE", result["blocker"])

    def test_evidence_refinement_preserves_bracket_on_unknown(self) -> None:
        # Checker-backed endpoints exercise the protocol, NOT one Figure 13 family.
        from make_it_stand.verification.faster_infeasibility import (
            verify_axis_aligned_yoyo_moment_capacity,
            verify_centered_axis_ellipsoid_witness,
        )
        from make_it_stand.verification.feasibility import FeasibilityReport
        fixtures = ROOT / "experiments/spin_it_reproduction/faster/fixtures"
        feasible = verify_centered_axis_ellipsoid_witness(fixtures / "cherchi-centered-axis-fixed-shell-v1/problem.json")
        infeasible = verify_axis_aligned_yoyo_moment_capacity(fixtures / "cherchi-rotated-box-fixed-shell-v2/problem.json")
        unknown = FeasibilityReport(reason="NO_CONTINUOUS_PROOF", input_hash="", system_hash="")
        runner = load_runner()
        result = runner.refine_feasibility_bracket(
            lambda d: feasible if d == 0.0 else infeasible if d == 1.0 else unknown,
            lower=0.0, upper=1.0, tolerance=1e-11,
        )
        self.assertEqual(result["status"], "UNKNOWN")
        self.assertEqual(result["lower"], 0.0)
        self.assertEqual(result["upper"], 1.0)
        self.assertIsNone(result["critical_distance"])
        self.assertEqual(result["observations"][-1]["evidence"]["reason"], "NO_CONTINUOUS_PROOF")
        self.assertEqual(result["family_monotonicity"], "NOT_VERIFIED")

    def test_evidence_refinement_rejects_solver_boolean_and_discrete_scope(self) -> None:
        from make_it_stand.verification.feasibility import FeasibilityReport
        runner = load_runner()
        with self.assertRaises(TypeError):
            runner.refine_feasibility_bracket(lambda d: True, lower=0.0, upper=1.0)
        discrete = FeasibilityReport(reason="SOLVER_DID_NOT_CONVERGE", input_hash="", system_hash="")
        result = runner.refine_feasibility_bracket(lambda d: discrete, lower=0.0, upper=1.0)
        self.assertEqual(result["status"], "UNKNOWN")
        self.assertEqual(result["reason"], "CONTINUOUS_ENDPOINT_PROOFS_REQUIRED")

    def test_final_payload_has_timing_independent_scientific_hash(self) -> None:
        runner = load_runner()
        first = runner.finalize_payload({"schema": "x", "elapsed_seconds": 1.0})
        second = runner.finalize_payload({"schema": "x", "elapsed_seconds": 9.0})
        self.assertEqual(first["scientific_hash"], second["scientific_hash"])


class ParetoInputDispatchTests(unittest.TestCase):
    def test_versioned_preflight_reports_its_actual_snapshot_not_mutable_narrative(self):
        runner = load_runner()
        document = json.loads((runner.ROOT.parent/"sources/case_asset_readiness.json").read_text())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = b"Independent versioned input-gap evidence.\n"
            digest = hashlib.sha256(source).hexdigest()
            snapshot = root/f"snapshots/{digest}.md"
            snapshot.parent.mkdir(); snapshot.write_bytes(source)
            document.update(schema="spin-it-case-asset-readiness-v2", schema_version=2,
                source_audit_path=f"snapshots/{digest}.md", source_audit_sha256=digest)
            path = root/"readiness.json"; path.write_text(json.dumps(document))
            result = runner.formal_pareto_preflight(path)
            self.assertEqual(result["readiness_validation"]["status"], "VALIDATED_INVENTORY")
            self.assertEqual(result["source_audit_sha256"], digest)
            self.assertEqual(result["source_audit_binding"]["path"], str(snapshot.resolve()))
            self.assertFalse(result["solver_executed"])
            self.assertEqual(result["status"], "BLOCKED")
            snapshot.unlink()  # Private temporary fixture only.
            rejected = runner.formal_pareto_preflight(path)
            self.assertEqual(rejected["readiness_validation"]["exception_type"], "FileNotFoundError")
            self.assertIsNone(rejected["source_audit_sha256"])
            self.assertTrue(all(row["asset_status"] == "UNVERIFIED" for row in rejected["figures"]))

    def test_unreadable_inventory_returns_serializable_blocker(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            for path, error_type in ((Path(directory)/"missing.json", "FileNotFoundError"),
                                     (Path(directory), "IsADirectoryError")):
                with self.subTest(error_type=error_type):
                    try:
                        result = runner.formal_pareto_preflight(path)
                    except OSError as error:
                        self.fail(f"input failure lost persistent blocker: {error}")
                    saved = json.loads(json.dumps(result, allow_nan=False))
                    self.assertEqual(saved["status"], "BLOCKED")
                    self.assertFalse(saved["solver_executed"])
                    self.assertEqual(saved["readiness_validation"]["exception_type"], error_type)
                    self.assertIsNone(saved["readiness_sha256"])
                    self.assertEqual(saved["input_read_errors"][0]["exception_type"], error_type)
                    self.assertTrue(all(row["asset_status"] == "UNVERIFIED"
                                        for row in saved["figures"]))

    def test_cli_rejects_conflicting_paths_and_missing_static_weights(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/"out.json"
            choices = (["--run-pareto"], ["--pareto-static-case", "unused.json"],
                       ["--pareto-top-case", "unused.json"])
            arguments = [choices[i]+choices[j] for i,j in ((0,1),(0,2),(1,2))]
            arguments.append(choices[1])
            for args in arguments:
                with self.subTest(args=args):
                    result = subprocess.run([sys.executable, str(RUNNER), *args,
                        "--output", str(output)], capture_output=True, text=True, timeout=30)
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertFalse(output.exists())
                    self.assertNotIn("Traceback", result.stderr)

    def test_default_envelope_does_not_advertise_missing_originals_as_ready(self):
        payload = load_runner().build_payload()
        for key in ("figure11", "figure12"):
            self.assertEqual(payload[key]["status"], "BLOCKED_MISSING_ORIGINAL_FIGURE_INPUT")
            self.assertNotIn("result", payload[key])

    def test_original_pareto_cli_preserves_input_blocker_without_substitution(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/"report.json"
            result = subprocess.run([sys.executable,str(RUNNER),"--run-pareto",
                "--output",str(output)],capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,2,result.stderr)
            self.assertTrue(output.is_file(), "blocked original inputs need a persistent report")
            payload=json.loads(output.read_text())
            preflight=payload["formal_pareto_preflight"]
            self.assertEqual(preflight["status"],"BLOCKED")
            self.assertFalse(preflight["solver_executed"])
            rows={r["requirement_id"]:r for r in preflight["figures"]}
            self.assertEqual(set(rows),{"faster.figure.11","faster.figure.12"})
            self.assertTrue(all("blocker" in row for row in rows.values()))
            self.assertNotIn("result",payload["figure11"])
            self.assertNotIn("result",payload["figure12"])

    def test_stale_inventory_is_rejected_without_fabricating_missing_or_ready_assets(self):
        runner=load_runner()
        document=json.loads((runner.ROOT.parent/"sources/case_asset_readiness.json").read_text())
        document.update(schema="spin-it-case-asset-readiness-v1", schema_version=1)
        document.pop("source_audit_path", None)
        document["source_audit_sha256"]="0"*64
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"readiness.json";path.write_text(json.dumps(document))
            result=runner.formal_pareto_preflight(path)
            self.assertEqual(result["readiness_validation"]["status"],"REJECTED")
            self.assertIn("SHA-256 mismatch",result["readiness_validation"]["diagnostic"])
            self.assertTrue(all(row["asset_status"]=="UNVERIFIED" for row in result["figures"]))
            self.assertFalse(result["solver_executed"])

    def test_valid_inventory_retains_egg_gaps_and_cannot_itself_admit_a_request(self):
        runner=load_runner()
        document=json.loads((runner.ROOT.parent/"sources/case_asset_readiness.json").read_text())
        document.update(schema="spin-it-case-asset-readiness-v1", schema_version=1)
        document.pop("source_audit_path", None)
        # Isolated validator fixture, NOT a repair of the real inventory.
        document["source_audit_sha256"]=hashlib.sha256(
            (runner.ROOT.parent/"sources/strict_asset_audit.md").read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"readiness.json";path.write_text(json.dumps(document))
            result=runner.formal_pareto_preflight(path)
            self.assertEqual(result["readiness_validation"]["status"],"VALIDATED_INVENTORY")
            self.assertIn("egg"," ".join(result["figures"][0]["missing_artifacts"]))
            for row in document["cases"]:
                if row["requirement_id"] in ("faster.figure.11","faster.figure.12"):
                    row.update(asset_status="AUTHOR_ORIGINAL",missing_artifacts=[],missing_critical_parameters=[])
            path.write_text(json.dumps(document))
            result=runner.formal_pareto_preflight(path)
            self.assertEqual(result["status"],"BLOCKED")
            self.assertFalse(result["source_authenticity_verified"])
            self.assertFalse(result["solver_executed"])
            self.assertTrue(all(row["blocker"]=="MISSING_AUTHENTICATED_LOCKED_FIGURE_REQUEST_BINDING"
                                for row in result["figures"]))

    def test_actual_static_method_cli_is_separate_and_repeatable(self):
        from test_faster_static_table_dispatch import method_case
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            case=method_case(root, maximum_iterations=0)
            path=root/"case.json";path.write_text(json.dumps(case))
            reports=[]
            for repeat in range(2):
                output=root/f"method-{repeat}.json"
                result=subprocess.run([sys.executable,str(RUNNER),
                    "--pareto-static-case",str(path),"--pareto-weights","0","0.5","1",
                    "--output",str(output)],capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stderr)
                payload=json.loads(output.read_text()); reports.append(payload)
                method=payload["static_pareto_method"]
                self.assertEqual(method["solver_family"],"STATIC_PLANE_FAMILY_V1")
                self.assertEqual(method["formal_paper_admission"],"NOT_GRANTED_BY_METHOD_RUNNER")
                self.assertEqual([r["weight"] for r in method["rows"]],[0.,.5,1.])
                self.assertEqual([r["application"]["static_problem"]["static_objective"]
                    for r in method["rows"]],
                    ["LOWER_CENTER_OF_MASS","PARETO_TRADEOFF","MAXIMUM_MASS"])
                self.assertGreater(method["cost_ledger"]["function_evaluations_all_solves"],0)
                self.assertNotIn("result",payload["figure11"])
                self.assertEqual(payload["figure11"]["status"],"BLOCKED_MISSING_ORIGINAL_FIGURE_INPUT")
            self.assertEqual(reports[0]["scientific_hash"],reports[1]["scientific_hash"])

    def test_static_method_does_not_silently_dispatch_a_top(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=root/"case.json"
            case.write_text(json.dumps({"id":"wrong-mode","mode":"SPINNING_TOP"}))
            result=subprocess.run([sys.executable,str(RUNNER),"--pareto-static-case",str(case),
                "--pareto-weights","0","1","--output",str(root/"out.json")],
                capture_output=True,text=True,timeout=30)
            self.assertNotEqual(result.returncode,0)
            self.assertIn("static method case requires STATIC_BALANCE",result.stderr)
            self.assertFalse((root/"out.json").exists())


class SpinItFasterPaperFigureNativeTests(unittest.TestCase):
    def test_invalid_quadric_side_witness_is_classified_in_solver_event(self) -> None:
        runner = load_runner()
        manifest = json.loads(
            (runner.ROOT / "paper_cases.json").read_text(encoding="utf-8")
        )
        case = next(row for row in manifest["cases"] if row["id"] == "heart-spinning-top")
        table = runner._table_module()
        request = table.build_request(case)
        design = table.read_obj(runner.ROOT / case["design_domain"]["path"])
        target_points = table.target_domain_points(case, design[0])
        seed = 2024 + 7919 * sum(map(ord, case["id"]))
        initial = table.generate_initial_quadrics(case["mode"], seed, 8, target_points)[0]

        from make_it_stand.spin import solve_spin_it_faster

        result = solve_spin_it_faster(request, initial_quadrics=(initial,))
        rejected = [
            event for event in result.events
            if event.reason == "evaluation_exception"
            and "outward quadric interface faces disagree on the interior sign"
            in event.detail
        ]
        self.assertTrue(rejected, "fixed source case must preserve its failed geometry trial")
        self.assertIn("positive_exterior_faces=", rejected[0].detail)
        self.assertIn("negative_exterior_faces=", rejected[0].detail)
        self.assertEqual(result.status, "CONVERGED")
        self.assertTrue(result.succeeded)
        self.assertTrue(result.geometry_sealed)

    def test_failed_pareto_reference_starts_remain_in_full_native_ledger(self):
        runner = load_runner()
        table = runner._table_module()
        cases = json.loads(table.DEFAULT_CASES.read_text())['cases']
        case = next(row for row in cases if row['id']=='static-tilted-vase')
        valid = table.generate_initial_quadrics(case['mode'], 0, 1,
            table.target_domain_points(case, table.read_case_mesh(case,'design_domain')[0]))
        invalid = [-1.] + [0.]*9
        # Freeze a deliberately invalid method-regression start, never use
        # this sequence as a paper protocol or replace ordinary failed seeds.
        with patch.object(runner, '_table_module', return_value=table), \
             patch.object(table, 'generate_initial_quadrics', return_value=[invalid, *valid]):
            result = runner.run_legacy_quadric_pareto_sweep(case, [0., .5, 1.])
        references = result['reference_attempts']
        self.assertEqual(references[0]['solve_record']['initial_quadric'], invalid)
        self.assertNotEqual(references[0]['solve_record']['status'], 'CONVERGED')
        self.assertEqual(references[0]['solve_record']['function_evaluations'], 1)
        records = [row['solve_record'] for row in references]
        records += [row['solve_record'] for row in result['rows']]
        self.assertEqual(result['cost_ledger']['function_evaluations_all_solves'],
                         sum(record['function_evaluations'] for record in records))
        for record in records:
            table.validate_attempt_ledger([record], expected_runs=1, maximum_restarts=1)
        if result['reference_status']!='CONVERGED':
            self.assertEqual(len(references), 2)
            self.assertEqual(result['rows'], [])
            self.assertFalse(result['front_audit']['finite'])

    def test_all_failed_pareto_reference_search_preserves_every_start(self):
        runner = load_runner()
        table = runner._table_module()
        case = next(row for row in json.loads(table.DEFAULT_CASES.read_text())['cases']
                    if row['id']=='static-tilted-vase')
        invalid = [-1.] + [0.]*9
        with patch.object(runner, '_table_module', return_value=table), \
             patch.object(table, 'generate_initial_quadrics', return_value=[invalid]*3):
            result = runner.run_legacy_quadric_pareto_sweep(case, [0., .5, 1.])
        self.assertEqual(len(result['reference_attempts']), 3)
        self.assertEqual(result['cost_ledger']['function_evaluations_all_solves'], 3)
        self.assertEqual(result['rows'], [])
        self.assertFalse(result['front_audit']['finite'])
        for row in result['reference_attempts']:
            self.assertNotEqual(row['solve_record']['status'], 'CONVERGED')
            self.assertFalse(row['solve_record']['numerical_acceptance']['passed'])
        # Explicit tags preserve nonfinite diagnostics without JSON NaN/Inf.
        json.dumps(result, allow_nan=False)

    def test_three_point_static_and_top_pareto_sweeps(self) -> None:
        runner = load_runner()
        manifest = json.loads(
            (ROOT / "experiments/spin_it_reproduction/faster/paper_cases.json").read_text(
                encoding="utf-8"
            )
        )
        cases = {case["id"]: case for case in manifest["cases"]}
        static = runner.run_legacy_quadric_pareto_sweep(
            cases["static-tilted-vase"],
            [0.0, 0.5, 1.0],
        )
        top = runner.run_legacy_quadric_pareto_sweep(
            cases["heart-spinning-top"],
            [0.0, 0.5, 1.0],
        )
        for sweep in (static, top):
            self.assertIn('reference_attempts', sweep)
            references = sweep['reference_attempts']
            self.assertTrue(references)
            self.assertEqual([row['candidate_index'] for row in references], list(range(len(references))))
            records = [row['solve_record'] for row in references]
            records += [row['solve_record'] for row in sweep['rows']]
            for record in records:
                self.assertEqual(sum(event['type']=='evaluation' for event in record['solver_events']),
                                 record['function_evaluations'])
                self.assertEqual(record['solver_events'][0]['type'], 'attempt_start')
                self.assertEqual(record['solver_events'][-1]['type'], 'attempt_terminal')
                # Preserve the supplied start and the native renormalized
                # start separately; ordinary unit normalization may round.
                actual_start = record['solver_events'][0]['quadric_coefficients']
                self.assertEqual(len(actual_start), len(record['initial_quadric']))
                self.assertLessEqual(max(abs(a-b) for a,b in
                    zip(record['initial_quadric'], actual_start)), 1e-14)
                self.assertIn('raw_equation_residual', record)
                self.assertIn('numerical_acceptance', record)
            self.assertEqual(sweep['cost_ledger']['function_evaluations_all_solves'],
                sum(row['function_evaluations'] for row in records))
            self.assertEqual(sweep['cost_ledger']['reference_attempts'], len(references))
            self.assertEqual(sweep['formal_paper_admission'], 'NOT_GRANTED_LEGACY_RECONSTRUCTED_DISPATCH')
            self.assertEqual(sweep["reference_status"], "CONVERGED")
            self.assertTrue(all(row["status"] == "CONVERGED" for row in sweep["rows"]))
            self.assertTrue(sweep["front_audit"]["finite"])
            self.assertTrue(sweep["front_audit"]["nondominated"])
        self.assertEqual(
            {row["interface_subtype"] for row in static["rows"]},
            {"PLANE"},
        )
        self.assertGreaterEqual(
            len({row["interface_subtype"] for row in top["rows"]}),
            2,
        )

    def test_tilted_cube_missing_proof_blocks_legacy_status_search(self) -> None:
        runner = load_runner()
        result = runner.run_tilted_cube_feasibility(tolerance=0.05, upper=0.8)
        self.assertEqual(result["asset_class"], "RECONSTRUCTED_ASSET")
        self.assertEqual(result["feasibility_status"], "UNKNOWN")
        self.assertIsNone(result["critical_distance"])
        self.assertEqual(result["rows"], [])


if __name__ == "__main__":
    unittest.main()
