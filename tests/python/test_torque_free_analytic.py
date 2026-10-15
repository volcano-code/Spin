"""Analytic torque-free oracles, not an original-paper experiment.

Breaks caught: dropping off-diagonal inertia, reversing quaternion composition,
wrong body/world convention, missing/invalid angular-velocity trace, and fake
conservation summaries. Expected trajectories use a closed form, not the ODE.
"""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.spatial.transform import Rotation


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'torque_free_subject', ROOT/'experiments/reconstructed_heart_top/heart_top.py')
SUBJECT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUBJECT)


def fixture():
    # Rotated axisymmetric inertia with all three off-diagonal terms nonzero.
    basis = Rotation.from_rotvec([.41, -.29, .23]).as_matrix()
    principal = np.array([.2, .2, .3])
    omega = np.array([10., 14., 40.*np.pi])
    return basis, principal, omega


def analytic(times, basis, principal, omega, azimuth, tilt=2.):
    """R(t)=Rot(L*t/I_perp) R(0) B Rot(e3*psi*t) B^T.

    In the principal frame Euler's equations give a transverse angular
    velocity rotating at k=(I_axial/I_perp-1)*omega_axial. The opposing
    right rotation psi=-k leaves spatial angular momentum fixed. This closed
    form needs two equal principal moments; it is not a triaxial solver.
    """
    if principal[0] != principal[1]:
        raise ValueError('axisymmetric reference only')
    times = np.asarray(times)
    initial = Rotation.from_rotvec(np.deg2rad(tilt)*np.array([
        np.cos(azimuth), np.sin(azimuth), 0.])).as_matrix()
    spatial_L = initial @ basis @ (principal*omega)
    k = (principal[2]/principal[0]-1.)*omega[2]
    left = Rotation.from_rotvec(times[:, None]*spatial_L/principal[0]).as_matrix()
    right = Rotation.from_rotvec(np.outer(times, [0., 0., -k])).as_matrix()
    rotations = left @ initial @ basis @ right @ basis.T
    cos, sin = np.cos(k*times), np.sin(k*times)
    wp = np.column_stack((omega[0]*cos-omega[1]*sin,
                          omega[0]*sin+omega[1]*cos,
                          np.full(len(times), omega[2])))
    return rotations, wp @ basis.T


def measure(result, inertia, expected_R, expected_w, *, expected_times):
    """Independent full-trace checks; ignore reported conservation booleans."""
    t = np.asarray(result['time'])
    grid = np.asarray(expected_times)
    if (grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all()
            or not np.all(np.diff(grid) > 0) or t.shape != grid.shape
            or not np.allclose(t, grid, rtol=0.,
                atol=32*np.finfo(float).eps*max(1., float(np.max(np.abs(grid)))))):
        raise ValueError('trajectory does not cover the independently requested time grid')
    q = np.asarray(result['quaternions'])
    w = np.asarray(result['angular_velocity_body'])
    if (t.ndim != 1 or len(t) < 2 or q.shape != (len(t), 4)
            or w.shape != (len(t), 3) or expected_R.shape != (len(t), 3, 3)
            or expected_w.shape != w.shape or not np.all(np.diff(t) > 0)
            or not all(np.isfinite(a).all() for a in (t, q, w, inertia, expected_R, expected_w))):
        raise ValueError('invalid or incomplete dynamics trace')
    norm_error = float(np.max(np.abs(np.linalg.norm(q, axis=1)-1)))
    R = Rotation.from_quat(q).as_matrix()
    orientation_error = float(np.max(Rotation.from_matrix(
        R @ np.swapaxes(expected_R, 1, 2)).magnitude()))
    omega_error = float(np.max(np.linalg.norm(w-expected_w, axis=1)) /
                        np.linalg.norm(expected_w[0]))
    energy = np.einsum('ni,ij,nj->n', w, inertia, w)/2
    momentum = np.einsum('nij,jk,nk->ni', R, inertia, w)
    return {
        'maximum_orientation_error_radians': orientation_error,
        'maximum_relative_omega_error': omega_error,
        'maximum_quaternion_norm_error': norm_error,
        'independent_relative_energy_drift': float(np.max(np.abs(energy-energy[0]))/energy[0]),
        'independent_relative_world_momentum_drift': float(np.max(
            np.linalg.norm(momentum-momentum[0], axis=1))/np.linalg.norm(momentum[0])),
    }


class TorqueFreeAnalyticTests(unittest.TestCase):
    def test_full_tensor_trajectory_matches_independent_closed_form(self):
        basis, principal, omega = fixture()
        I = basis @ np.diag(principal) @ basis.T
        self.assertTrue(np.all(np.abs(I[np.triu_indices(3, 1)]) > .001))
        result = SUBJECT.simulate(1., np.zeros(3), I, basis@omega,
                                  gravity=0., duration=.125, azimuth=.7)
        self.assertIn('angular_velocity_body', result,
                      'independent invariants require the actual integrated omega trace')
        R, w = analytic(result['time'], basis, principal, omega, .7)
        for name, value in measure(result, I, R, w, expected_times=np.arange(31)/240).items():
            self.assertLess(value, 1e-8, name)

    def test_spherical_reference_reduces_to_constant_body_spin(self):
        times = np.array([0., .13, .71])
        B = Rotation.from_rotvec([.7, .2, -.4]).as_matrix()
        wp = np.array([2., -1., .3])
        R, w = analytic(times, B, np.ones(3), wp, 0., tilt=0.)
        expected = Rotation.from_rotvec(times[:, None]*(B@wp)).as_matrix()
        np.testing.assert_allclose(R, expected, atol=1e-14)
        np.testing.assert_allclose(w, np.tile(B@wp, (len(times), 1)), atol=1e-14)

    def test_checker_detects_dropped_cross_inertia(self):
        B, principal, wp = fixture()
        I = B@np.diag(principal)@B.T
        wrong = SUBJECT.simulate(1., np.zeros(3), np.diag(np.diag(I)), B@wp,
                                 gravity=0., duration=.125, azimuth=.7)
        self.assertIn('angular_velocity_body', wrong)
        R, w = analytic(wrong['time'], B, principal, wp, .7)
        self.assertGreater(measure(wrong, I, R, w, expected_times=np.arange(31)/240)
                           ['maximum_orientation_error_radians'], .01)

    def test_checker_rejects_truncation_and_nan(self):
        R = np.tile(np.eye(3), (2, 1, 1)); w = np.ones((2, 3))
        valid = dict(time=[0., 1.], quaternions=[[0., 0., 0., 1.]]*2,
                     angular_velocity_body=w.tolist())
        with self.assertRaises(ValueError):
            measure(dict(valid, angular_velocity_body=[[1., 1., 1.]]), np.eye(3), R, w,
                    expected_times=np.array([0., 1.]))
        with self.assertRaises(ValueError):
            measure(dict(valid, quaternions=[[0., 0., 0., 1.], [0., 0., 0., float('nan')]]),
                    np.eye(3), R, w, expected_times=np.array([0., 1.]))

    def test_checker_rejects_consistently_truncated_and_irregular_time_grids(self):
        for times in ([0., .5], [0., .6, 1.]):
            with self.subTest(times=times):
                n = len(times)
                R = np.tile(np.eye(3), (n, 1, 1)); w = np.ones((n, 3))
                result = dict(time=times, quaternions=[[0., 0., 0., 1.]]*n,
                              angular_velocity_body=w.tolist())
                with self.assertRaises(ValueError):
                    measure(result, np.eye(3), R, w, expected_times=np.array([0., .5, 1.]))


if __name__ == '__main__':
    unittest.main()
