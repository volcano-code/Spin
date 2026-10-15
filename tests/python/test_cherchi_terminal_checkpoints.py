"""A terminal audit interruption must not erase the actual accepted material."""
import hashlib
import json
from dataclasses import replace

import pytest

from make_it_stand import spin
from make_it_stand.verification import cherchi_inprocess as ci
from make_it_stand.verification import faster_material as fm
from test_cherchi_rational_boundary_profile import CORE, actual_material


def settings():
    config = json.loads((CORE/'core-benchmark-inputs-job-01/output/small.json').read_text())
    return spin.SpinEvaluatorConfig(**config['evaluator']), spin.SpinSolverConfig(**config['request']['settings'])


def terminal_solver(monkeypatch, q):
    # Replace only Newton iteration, retaining real input preparation, Cherchi,
    # material validation, journal publication and terminal control flow.
    def solve(request, evaluate, **kwargs):
        mapped = evaluate(q)
        return dict(status='CONVERGED', quadric_coefficients=q,
                    raw_equation_residual=1e-14, events=[],
                    trace=[dict(material_geometry_hash=mapped['geometry_hash'])])
    monkeypatch.setattr(fm, 'solve_spin_it_faster_with_external_moment_evaluator', solve)


def test_numerical_terminal_is_durable_before_complement_interrupt(tmp_path, monkeypatch, actual_material):
    problem, _, q, _ = actual_material
    terminal_solver(monkeypatch, q)
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt('interrupted after Newton')
    monkeypatch.setattr(ci.InProcessCherchiEvaluator, 'complement', interrupt)
    evaluator, solver = settings()
    output = tmp_path/'run'
    with pytest.raises(KeyboardInterrupt):
        ci.solve_cherchi_inprocess_and_seal_faster(problem,
            evaluator_settings=evaluator, solver_settings=solver, output_directory=output,
            rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
            solver_boundary_profile='EXACT_RATIONAL_BOUNDARY_V1')
    assert (output/'native-result.json').is_file(), 'Newton result was lost before sealing'
    terminal = json.loads((output/'terminal-material-reference.json').read_text())
    material_path = output/terminal['evaluation_file']
    assert terminal['coefficients'] == q
    assert terminal['sha256'] == hashlib.sha256(material_path.read_bytes()).hexdigest()
    assert not (output/'complete.json').exists()


def test_complement_is_durable_before_independent_audit_interrupt(tmp_path, monkeypatch, actual_material):
    problem, _, q, _ = actual_material
    terminal_solver(monkeypatch, q)
    original = fm._external_geometry_moment_map
    def interrupt_complement(*args, **kwargs):
        if kwargs.get('require_derivative') is False:
            raise KeyboardInterrupt('interrupted before complement audit')
        return original(*args, **kwargs)
    monkeypatch.setattr(fm, '_external_geometry_moment_map', interrupt_complement)
    evaluator, solver = settings()
    output = tmp_path/'run'
    with pytest.raises(KeyboardInterrupt):
        ci.solve_cherchi_inprocess_and_seal_faster(problem,
            evaluator_settings=evaluator, solver_settings=solver, output_directory=output,
            rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
            solver_boundary_profile='EXACT_RATIONAL_BOUNDARY_V1')
    assert (output/'terminal-complement.json').is_file(), 'same-call complement was lost before audit'
    complement = json.loads((output/'terminal-complement.json').read_text())
    retained = json.loads((output/'evaluation-000001.json').read_text())
    assert complement['quadric_source_identity'] == retained['quadric_source_identity']
    assert complement['boolean_operation'] != retained['boolean_operation']
    assert not (output/'complete.json').exists()


def test_formal_run_prepares_locked_input_only_once(tmp_path, monkeypatch, actual_material):
    problem, _, q, _ = actual_material
    terminal_solver(monkeypatch, q)
    calls = []
    original = fm._prepare
    def counted(*args, **kwargs):
        calls.append(args[0])
        return original(*args, **kwargs)
    monkeypatch.setattr(fm, '_prepare', counted)
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt('stop after observing preparation')
    monkeypatch.setattr(ci.InProcessCherchiEvaluator, 'complement', interrupt)
    evaluator, solver = settings()
    with pytest.raises(KeyboardInterrupt):
        ci.solve_cherchi_inprocess_and_seal_faster(problem,
            evaluator_settings=evaluator, solver_settings=solver, output_directory=tmp_path/'run',
            rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
            solver_boundary_profile='EXACT_RATIONAL_BOUNDARY_V1')
    assert len(calls) == 1, 'formal entry repeated input partition and fixed integration'


@pytest.mark.parametrize('fault', ['dictionary', 'different_seed', 'different_input'])
def test_prepared_problem_cannot_bypass_input_or_settings_validation(actual_material, fault):
    problem, _, q, material = actual_material
    evaluator, solver = settings()
    prepared = fm._prepare_faster_problem(problem, evaluator, solver)
    problem = prepared.problem
    if fault == 'dictionary':
        prepared = dict(problem=problem, passed=True)
    elif fault == 'different_seed':
        solver = replace(solver, deterministic_seed=solver.deterministic_seed+1)
    else:
        problem = actual_material[0]
    with pytest.raises(ValueError, match='prepared Faster problem'):
        fm.solve_cherchi_geometry_and_seal_faster(problem,
            geometry_evaluator=lambda *args: material, complement_evaluator=lambda *args: material,
            initial_quadrics=(), evaluator_settings=evaluator, solver_settings=solver,
            _prepared_problem=prepared)


def test_complement_recovery_uses_saved_source_without_newton_or_source_reconstruction(actual_material, monkeypatch):
    problem, request, q, material = actual_material
    evaluator_settings, _ = settings()
    producer = ci.InProcessCherchiEvaluator(request, evaluator_settings,
        normalization=problem.normalization,
        design_vertices_physical=problem.design_regions[0].mesh.arrays()[0],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    def forbidden(*args, **kwargs):
        raise AssertionError('recovery must not optimize or reconstruct its accepted source')
    monkeypatch.setattr(ci.InProcessCherchiEvaluator, '__call__', forbidden)
    monkeypatch.setattr(fm, 'solve_spin_it_faster_with_external_moment_evaluator', forbidden)
    recover = getattr(producer, 'complement_of_saved_source', None)
    assert callable(recover), 'missing saved-source-only complement recovery'
    packets = []
    complement = recover(material, packet_sink=packets.append)
    assert len(packets) == 1
    assert packets[0]['exact_vertices'] == complement['material_proof']['exact_vertices']
    assert producer.call_counts['source'] == 0
    assert producer.call_counts['boolean'] == 1
    assert complement['quadric_coefficients'] == tuple(q)
    assert complement['quadric_source_identity'] == material['quadric_source_identity']
