import numpy as np

from rpr_model import GAMMA, GAMMA5, I4, ModelParameters, RPRModel, dirac_spinor, minkowski_dot


def test_clifford_algebra():
    metric = [1.0, -1.0, -1.0, -1.0]
    for mu in range(4):
        for nu in range(4):
            anti = GAMMA[mu] @ GAMMA[nu] + GAMMA[nu] @ GAMMA[mu]
            expected = 2.0 * metric[mu] * I4 if mu == nu else np.zeros((4, 4), complex)
            np.testing.assert_allclose(anti, expected, atol=1e-13)
    np.testing.assert_allclose(GAMMA5 @ GAMMA5, I4, atol=1e-13)


def test_spinor_normalization_and_dirac_equation():
    mass = 1.32171
    p = np.array([np.sqrt(mass**2 + 0.4**2), 0.3, 0.0, -np.sqrt(0.4**2 - 0.3**2)])
    for spin in (0, 1):
        u = dirac_spinor(p, mass, spin)
        ubar_u = u.conj().T @ GAMMA[0] @ u
        np.testing.assert_allclose(ubar_u, 1.0, atol=1e-13)


def test_kinematic_invariants():
    model = RPRModel()
    kin = model.kinematics(2.15, 0.91)
    np.testing.assert_allclose(minkowski_dot(kin.k1, kin.k1), model.par.m_k**2, atol=1e-12)
    np.testing.assert_allclose(minkowski_dot(kin.p1, kin.p1), model.par.m_p**2, atol=1e-12)
    np.testing.assert_allclose(minkowski_dot(kin.k2, kin.k2), model.par.m_k**2, atol=1e-12)
    np.testing.assert_allclose(minkowski_dot(kin.p2, kin.p2), model.par.m_xi**2, atol=1e-12)
    np.testing.assert_allclose(kin.k1 + kin.p1, kin.k2 + kin.p2, atol=1e-12)


def test_polarization_is_bounded_and_collinear_zero():
    model = RPRModel(ModelParameters(resonance_form="gaussian_squared"))
    for cth in np.linspace(-1.0, 1.0, 21):
        obs = model.observables(2.15, float(cth))
        assert obs["dsigma_domega"] >= 0.0
        assert abs(obs["polarization_y"]) <= 1.0 + 1e-12
    assert abs(model.observables(2.15, 1.0)["polarization_y"]) < 1e-12
    assert abs(model.observables(2.15, -1.0)["polarization_y"]) < 1e-12


def test_component_sum():
    model = RPRModel()
    components = model.amplitude_components(2.15, 0.8)
    summed = sum(components.values(), np.zeros((2, 2), complex))
    np.testing.assert_allclose(summed, model.amplitude(2.15, 0.8), atol=1e-13)

