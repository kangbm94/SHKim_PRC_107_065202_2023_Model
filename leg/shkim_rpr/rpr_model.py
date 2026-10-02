"""Reference implementation of the u+s part of Kim et al., PRC 107, 065202.

The implementation follows Eqs. (6)--(16) of the paper for
K- p -> K+ Xi-.  Dirac spinors are normalized as ubar u = 1, matching the
statement below Eq. (8).  The model deliberately exposes conventions that are
not fully specified in the article (form-factor power and resonance Gaussian)
so that they can be selected by validation against the published curves.

This module does not yet include the meson-baryon rescattering loop of
Eqs. (17)--(28).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import pi, sqrt
from typing import Dict, Iterable, Mapping, Optional

import numpy as np
from numpy.typing import NDArray
from scipy.special import gamma as gamma_function


ComplexMatrix = NDArray[np.complex128]
FourVector = NDArray[np.float64]


GEV2_TO_MICROBARN = 389.3793721


def _dirac_matrices() -> tuple[list[ComplexMatrix], ComplexMatrix, ComplexMatrix]:
    """Return gamma^mu, gamma5, and the 4x4 identity in the Dirac basis."""

    i2 = np.eye(2, dtype=np.complex128)
    z2 = np.zeros((2, 2), dtype=np.complex128)
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)

    g0 = np.block([[i2, z2], [z2, -i2]])
    gamma = [g0]
    for sigma in (sx, sy, sz):
        gamma.append(np.block([[z2, sigma], [-sigma, z2]]))
    gamma5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
    return gamma, gamma5, np.eye(4, dtype=np.complex128)


GAMMA, GAMMA5, I4 = _dirac_matrices()
METRIC = np.array([1.0, -1.0, -1.0, -1.0])
GAMMA_LOWER = [METRIC[mu] * GAMMA[mu] for mu in range(4)]
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)


def minkowski_dot(a: FourVector, b: FourVector) -> float:
    return float(np.sum(METRIC * a * b))


def lower(a: FourVector) -> FourVector:
    return METRIC * a


def slash(a: FourVector) -> ComplexMatrix:
    a_lower = lower(a)
    return sum((GAMMA[mu] * a_lower[mu] for mu in range(4)), np.zeros((4, 4), complex))


def dirac_spinor(momentum: FourVector, mass: float, spin: int) -> NDArray[np.complex128]:
    """Positive-energy spinor with ubar u = 1 and z-quantized spin."""

    if spin not in (0, 1):
        raise ValueError("spin must be 0 (up) or 1 (down)")
    energy = float(momentum[0])
    p = momentum[1:]
    chi = np.array([1.0, 0.0], complex) if spin == 0 else np.array([0.0, 1.0], complex)
    sigma_p = np.array(
        [[p[2], p[0] - 1j * p[1]], [p[0] + 1j * p[1], -p[2]]],
        dtype=np.complex128,
    )
    norm = sqrt((energy + mass) / (2.0 * mass))
    return norm * np.concatenate((chi, sigma_p @ chi / (energy + mass)))


def kallen(x: float, y: float, z: float) -> float:
    return x * x + y * y + z * z - 2.0 * (x * y + y * z + z * x)


@dataclass(frozen=True)
class Kinematics:
    sqrt_s: float
    cos_theta: float
    k1: FourVector
    p1: FourVector
    k2: FourVector
    p2: FourVector
    q_s: FourVector
    q_u: FourVector
    s: float
    u: float
    p_initial: float
    p_final: float


@dataclass(frozen=True)
class Resonance:
    name: str
    family: str
    spin: float
    parity: int
    mass: float
    width: float
    g_kn: float
    g_kxi: float


@dataclass
class ModelParameters:
    # Masses in GeV.
    m_k: float = 0.493677
    m_p: float = 0.9382720813
    m_xi: float = 1.32171
    m_lambda: float = 1.115683
    m_sigma: float = 1.192642
    m_sigma1385: float = 1.3837

    # Ground-state SU(3) couplings, Eq. (10).
    g_kn_lambda: float = -13.24
    g_kxi_lambda: float = 3.52
    g_kn_sigma: float = 3.58
    g_kxi_sigma: float = -13.26
    g_kn_sigma1385: float = -3.22
    g_kxi_sigma1385: float = -3.22

    # Regge scale factors, Eqs. (6)--(7).
    eta_lambda: float = 2.6
    eta_sigma: float = 0.66
    eta_sigma1385: float = 0.66
    regge_cutoff: float = 1.0
    regge_s0: float = 1.0

    # s-channel form factors, Eqs. (13) and (16).
    born_cutoff: float = 0.85
    born_n: int = 2
    resonance_cutoff: float = 0.73
    # Validation against the isolated curves in Fig. 4 shows that the quoted
    # F_Y is the diagram-level factor used once in the published calculation.
    form_factor_power: int = 1
    # "paper_literal": exp[-(s-M^2)/Lambda^2]
    # "gaussian_squared": exp[-(s-M^2)^2/Lambda^4]
    resonance_form: str = "gaussian_squared"
    # High-spin s-channel projector convention.  Ref. 46 uses sqrt(s) for
    # above-threshold resonances; "pole_mass" keeps the literal mass argument.
    s_projector_mass_mode: str = "pole_mass"

    # Charged-channel isospin factors.  The relative Lambda/Sigma sign in the
    # s channel follows the field conventions in Eqs. (1)--(2).
    isospin_s_lambda: float = -1.0
    isospin_s_sigma: float = 1.0
    isospin_u_lambda: float = -1.0
    isospin_u_sigma: float = 1.0
    isospin_u_sigma1385: float = 1.0

    # Multiplicative amplitude scales used by the fitter.  A value of one is
    # the published parameter set.  The overall normalization is separate so
    # polarization-only fits cannot spuriously determine it.
    scales: Dict[str, float] = field(default_factory=lambda: {
        "u_lambda": 1.0,
        "u_sigma": 1.0,
        "u_sigma1385": 1.0,
        "s_lambda": 1.0,
        "s_sigma": 1.0,
        "Lambda1890": 1.0,
        "Lambda2100": 1.0,
        "Sigma2030": 1.0,
        "Sigma2230": 1.0,
        "Sigma2250": 1.0,
    })
    cross_section_normalization: float = 1.0

    resonances: tuple[Resonance, ...] = (
        Resonance("Lambda1890", "Lambda", 1.5, +1, 1.890, 0.120, 0.84, -0.26),
        Resonance("Lambda2100", "Lambda", 3.5, -1, 2.100, 0.200, 2.41, 2.95),
        Resonance("Sigma2030", "Sigma", 3.5, +1, 2.030, 0.180, 0.82, -0.93),
        Resonance("Sigma2230", "Sigma", 1.5, +1, 2.230, 0.345, 0.41, 0.34),
        # The article sets the mass to 2.290 GeV, despite the Sigma(2250) label.
        Resonance("Sigma2250", "Sigma", 3.5, -1, 2.290, 0.100, 0.59, 0.88),
    )


class RPRModel:
    """Free-proton K- p -> K+ Xi- RPR calculation without rescattering."""

    def __init__(self, parameters: Optional[ModelParameters] = None):
        self.par = parameters or ModelParameters()

    def kinematics(self, sqrt_s: float, cos_theta: float) -> Kinematics:
        p = self.par
        s = sqrt_s * sqrt_s
        if sqrt_s < p.m_k + p.m_xi:
            raise ValueError("sqrt(s) is below the K Xi threshold")
        if not -1.0 <= cos_theta <= 1.0:
            raise ValueError("cos(theta) must lie in [-1,1]")

        pi2 = kallen(s, p.m_k**2, p.m_p**2) / (4.0 * s)
        pf2 = kallen(s, p.m_k**2, p.m_xi**2) / (4.0 * s)
        p_i, p_f = sqrt(max(0.0, pi2)), sqrt(max(0.0, pf2))
        sin_theta = sqrt(max(0.0, 1.0 - cos_theta * cos_theta))

        e_ki = sqrt(p.m_k**2 + p_i**2)
        e_pi = sqrt(p.m_p**2 + p_i**2)
        e_kf = sqrt(p.m_k**2 + p_f**2)
        e_xf = sqrt(p.m_xi**2 + p_f**2)

        k1 = np.array([e_ki, 0.0, 0.0, p_i])
        p1 = np.array([e_pi, 0.0, 0.0, -p_i])
        k2 = np.array([e_kf, p_f * sin_theta, 0.0, p_f * cos_theta])
        p2 = np.array([e_xf, -p_f * sin_theta, 0.0, -p_f * cos_theta])
        q_s = k1 + p1
        q_u = p2 - k1
        u = minkowski_dot(q_u, q_u)
        return Kinematics(
            sqrt_s, cos_theta, k1, p1, k2, p2, q_s, q_u, s, u, p_i, p_f
        )

    @staticmethod
    def _contract_p(a: FourVector, b: FourVector, q: FourVector, mass: float) -> float:
        ql = lower(q)
        result = 0.0
        for mu in range(4):
            for nu in range(4):
                p_mn = -(METRIC[mu] if mu == nu else 0.0) + ql[mu] * ql[nu] / mass**2
                result += a[mu] * p_mn * b[nu]
        return float(result)

    @staticmethod
    def _contract_r(
        a: FourVector, b: FourVector, q: FourVector, mass: float
    ) -> ComplexMatrix:
        ql = lower(q)
        result = np.zeros((4, 4), dtype=np.complex128)
        for mu in range(4):
            for nu in range(4):
                r_mn = (
                    -GAMMA_LOWER[mu] @ GAMMA_LOWER[nu]
                    - (GAMMA_LOWER[mu] * ql[nu] - GAMMA_LOWER[nu] * ql[mu]) / mass
                    + I4 * ql[mu] * ql[nu] / mass**2
                )
                result += a[mu] * b[nu] * r_mn
        return result

    def _contract_delta3(
        self, a: FourVector, b: FourVector, q: FourVector, mass: float
    ) -> ComplexMatrix:
        # Eq. (9) is equivalently P_{mu nu} - R_{mu nu}/3.
        return self._contract_p(a, b, q, mass) * I4 - self._contract_r(a, b, q, mass) / 3.0

    def _contract_delta7(
        self, a: FourVector, b: FourVector, q: FourVector, mass: float
    ) -> ComplexMatrix:
        """Momentum-contracted spin-7/2 Behrends-Fronsdal projector.

        Because the three momenta at either vertex are identical, the 6- and
        36-term permutation sums reduce to four scalar/matrix contractions.
        This is the contracted form of the projector cited in Refs. 43--46.
        """

        a_pb = self._contract_p(a, b, q, mass)
        a_pa = self._contract_p(a, a, q, mass)
        b_pb = self._contract_p(b, b, q, mass)
        a_rb = self._contract_r(a, b, q, mass)
        scalar = a_pb**3 - (3.0 / 7.0) * a_pa * b_pb * a_pb
        matrix_coefficient = -(3.0 / 7.0) * a_pb**2 + (3.0 / 35.0) * a_pa * b_pb
        return scalar * I4 + matrix_coefficient * a_rb

    @staticmethod
    def _spin_matrix(operator: ComplexMatrix, kin: Kinematics, m_p: float, m_xi: float) -> ComplexMatrix:
        result = np.zeros((2, 2), dtype=np.complex128)
        for sf in range(2):
            uf = dirac_spinor(kin.p2, m_xi, sf)
            ubar = uf.conj().T @ GAMMA[0]
            for si in range(2):
                ui = dirac_spinor(kin.p1, m_p, si)
                result[sf, si] = ubar @ operator @ ui
        return result

    def _born_form_factor(self, s: float, mass: float) -> float:
        p = self.par
        numerator = p.born_n * p.born_cutoff**4
        base = numerator / (numerator + (s - mass**2) ** 2)
        return float(base**p.born_n)

    def _resonance_form_factor(self, s: float, mass: float) -> float:
        p = self.par
        x = (s - mass**2) / p.resonance_cutoff**2
        if p.resonance_form == "paper_literal":
            return float(np.exp(-x))
        if p.resonance_form == "gaussian_squared":
            return float(np.exp(-(x * x)))
        raise ValueError(f"unknown resonance_form: {p.resonance_form}")

    def _regge_factor(
        self, s: float, u: float, intercept: float, slope: float, spin_offset: float, eta: float
    ) -> float:
        p = self.par
        alpha = intercept + slope * u
        # Eq. (7): the square applies to eta_H and to the monopole ratio.
        c_h = (eta * p.regge_cutoff**2 / (p.regge_cutoff**2 - u)) ** 2
        return float(
            c_h
            * (s / p.regge_s0) ** (alpha - spin_offset)
            * gamma_function(spin_offset - alpha)
            * slope
        )

    def amplitude_components(self, sqrt_s: float, cos_theta: float) -> Dict[str, ComplexMatrix]:
        p = self.par
        k = self.kinematics(sqrt_s, cos_theta)
        components: Dict[str, ComplexMatrix] = {}

        # Reggeized spin-1/2 Lambda and Sigma, Eqs. (6) and (8).
        for name, mass, g1, g2, iso, intercept, slope, eta in (
            ("u_lambda", p.m_lambda, p.g_kn_lambda, p.g_kxi_lambda, p.isospin_u_lambda, -0.65, 0.94, p.eta_lambda),
            ("u_sigma", p.m_sigma, p.g_kn_sigma, p.g_kxi_sigma, p.isospin_u_sigma, -0.79, 0.87, p.eta_sigma),
        ):
            prefactor = iso * g1 * g2 / ((p.m_p + mass) * (p.m_xi + mass))
            operator = slash(k.k1) @ GAMMA5 @ (slash(k.q_u) + mass * I4) @ slash(k.k2) @ GAMMA5
            regge = self._regge_factor(k.s, k.u, intercept, slope, 0.5, eta)
            components[name] = p.scales[name] * prefactor * regge * self._spin_matrix(operator, k, p.m_p, p.m_xi)

        # Reggeized Sigma(1385), Eqs. (6) and (8).
        delta3_u = self._contract_delta3(k.k1, k.k2, k.q_u, p.m_sigma1385)
        operator = (slash(k.q_u) + p.m_sigma1385 * I4) @ delta3_u
        prefactor = (
            p.isospin_u_sigma1385
            * p.g_kn_sigma1385
            * p.g_kxi_sigma1385
            / p.m_k**2
        )
        regge = self._regge_factor(k.s, k.u, -0.27, 0.9, 1.5, p.eta_sigma1385)
        components["u_sigma1385"] = (
            p.scales["u_sigma1385"]
            * prefactor
            * regge
            * self._spin_matrix(operator, k, p.m_p, p.m_xi)
        )

        # Ground-state s-channel Lambda and Sigma, Eqs. (12)--(13).
        for name, mass, g1, g2, iso in (
            ("s_lambda", p.m_lambda, p.g_kn_lambda, p.g_kxi_lambda, p.isospin_s_lambda),
            ("s_sigma", p.m_sigma, p.g_kn_sigma, p.g_kxi_sigma, p.isospin_s_sigma),
        ):
            prefactor = iso * g1 * g2 / (
                (k.s - mass**2) * (p.m_p + mass) * (p.m_xi + mass)
            )
            operator = slash(k.k2) @ GAMMA5 @ (slash(k.q_s) + mass * I4) @ slash(k.k1) @ GAMMA5
            ff = self._born_form_factor(k.s, mass) ** p.form_factor_power
            components[name] = p.scales[name] * prefactor * ff * self._spin_matrix(operator, k, p.m_p, p.m_xi)

        # Selected s-channel resonances, Eqs. (14)--(16) and Table II.
        for resonance in p.resonances:
            iso = p.isospin_s_lambda if resonance.family == "Lambda" else p.isospin_s_sigma
            projector_mass = k.sqrt_s if p.s_projector_mass_mode == "sqrt_s" else resonance.mass
            if resonance.spin == 1.5:
                delta = self._contract_delta3(k.k2, k.k1, k.q_s, projector_mass)
                vertex_gamma = I4 if resonance.parity > 0 else GAMMA5
                operator = (
                    vertex_gamma
                    @ (slash(k.q_s) + resonance.mass * I4)
                    @ delta
                    @ vertex_gamma
                    / p.m_k**2
                )
            elif resonance.spin == 3.5:
                delta = self._contract_delta7(k.k2, k.k1, k.q_s, projector_mass)
                vertex_gamma = I4 if resonance.parity > 0 else GAMMA5
                operator = (
                    vertex_gamma
                    @ (slash(k.q_s) + resonance.mass * I4)
                    @ delta
                    @ vertex_gamma
                    / p.m_k**6
                )
            else:
                raise NotImplementedError(f"spin {resonance.spin} is not selected in this paper")

            denominator = k.s - resonance.mass**2 + 1j * resonance.mass * resonance.width
            ff = self._resonance_form_factor(k.s, resonance.mass) ** p.form_factor_power
            prefactor = iso * resonance.g_kn * resonance.g_kxi / denominator
            components[resonance.name] = (
                p.scales[resonance.name]
                * prefactor
                * ff
                * self._spin_matrix(operator, k, p.m_p, p.m_xi)
            )

        return components

    def amplitude(
        self,
        sqrt_s: float,
        cos_theta: float,
        include: Optional[Iterable[str]] = None,
    ) -> ComplexMatrix:
        components = self.amplitude_components(sqrt_s, cos_theta)
        names = components.keys() if include is None else include
        return sum((components[name] for name in names), np.zeros((2, 2), complex))

    def observables(
        self,
        sqrt_s: float,
        cos_theta: float,
        include: Optional[Iterable[str]] = None,
    ) -> Mapping[str, float]:
        p = self.par
        kin = self.kinematics(sqrt_s, cos_theta)
        amp = self.amplitude(sqrt_s, cos_theta, include)
        density = amp @ amp.conj().T
        spin_sum = 0.5 * float(np.trace(density).real)
        prefactor = (
            p.m_p
            * p.m_xi
            / (16.0 * pi**2 * kin.s)
            * kin.p_final
            / kin.p_initial
            * GEV2_TO_MICROBARN
            * p.cross_section_normalization
        )
        dsigma = prefactor * spin_sum
        py = 0.0 if spin_sum == 0.0 else float(np.trace(density @ PAULI_Y).real / np.trace(density).real)
        return {
            "dsigma_domega": dsigma,
            "polarization_y": py,
            "py_dsigma_domega": py * dsigma,
        }

    def total_cross_section(
        self, sqrt_s: float, include: Optional[Iterable[str]] = None, order: int = 96
    ) -> float:
        nodes, weights = np.polynomial.legendre.leggauss(order)
        integral = sum(
            weight * self.observables(sqrt_s, float(node), include)["dsigma_domega"]
            for node, weight in zip(nodes, weights)
        )
        return float(2.0 * pi * integral)


def beam_momentum_to_sqrt_s(p_lab: float, m_k: float = 0.493677, m_p: float = 0.9382720813) -> float:
    e_k = sqrt(p_lab * p_lab + m_k * m_k)
    return sqrt(m_k * m_k + m_p * m_p + 2.0 * m_p * e_k)
