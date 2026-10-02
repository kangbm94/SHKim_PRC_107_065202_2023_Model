#!/usr/bin/env python3
"""Generate validation plots corresponding to Figs. 4(a), 5(a), and 7(a)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rpr_model import GEV2_TO_MICROBARN, PAULI_Y, RPRModel


U_CHANNEL = ("u_lambda", "u_sigma", "u_sigma1385")
#U_CHANNEL = ("u_lambda")
S_CHANNEL = (
    "s_lambda",
    "s_sigma",
    "Lambda1890",
    "Lambda2100",
    "Sigma2030",
    "Sigma2230"
    "Sigma2250",
)


def observables_from_components(model, kin, components, include=None):
    names = components.keys() if include is None else include
    amplitude = sum((components[name] for name in names), np.zeros((2, 2), complex))
    density = amplitude @ amplitude.conj().T
    trace = float(np.trace(density).real)
    p = model.par
    prefactor = (
        p.m_p * p.m_xi / (16.0 * np.pi**2 * kin.s)
        * kin.p_final / kin.p_initial * GEV2_TO_MICROBARN
    )
    dsigma = prefactor * 0.5 * trace
    py = 0.0 if trace == 0 else float(np.trace(density @ PAULI_Y).real / trace)
    return dsigma, py * dsigma


def figure4(model: RPRModel, output: Path) -> None:
    energies = np.linspace(model.par.m_k + model.par.m_xi + 1e-4, 3.05, 110)
    selections = {
        "u-channel Regge": U_CHANNEL,
        r"$\Lambda(1890)3/2^+$": ("Lambda1890",),
        r"$\Lambda(2100)7/2^-$": ("Lambda2100",),
        r"$\Sigma(2030)7/2^+$": ("Sigma2030",),
        r"$\Sigma(2230)3/2^+$": ("Sigma2230",),
        r"$\Sigma(2250)7/2^-$": ("Sigma2250",),
        "full u+s": None,
    }
    curves = {label: [] for label in selections}
    nodes, weights = np.polynomial.legendre.leggauss(24)
    for energy in energies:
        integrals = {label: 0.0 for label in selections}
        for cth, weight in zip(nodes, weights):
            kin = model.kinematics(float(energy), float(cth))
            components = model.amplitude_components(float(energy), float(cth))
            for label, include in selections.items():
                dsigma, _ = observables_from_components(model, kin, components, include)
                integrals[label] += weight * dsigma
        for label in selections:
            curves[label].append(2.0 * np.pi * integrals[label])
    styles = ["--", ":", "--", "--", ":", "-.", "-"]
    widths = [1.4, 1.2, 1.3, 1.3, 1.2, 1.2, 2.0]
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for (label, values), style, width in zip(curves.items(), styles, widths):
        ax.plot(energies, values, style, lw=width, label=label)
    ax.set(xlim=(1.8, 3.05), ylim=(0, 250), xlabel=r"$\sqrt{s}$ [GeV]", ylabel=r"$\sigma$ [$\mu$b]")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, fontsize=9, ncol=2)
    ax.set_title(r"$K^-p\rightarrow K^+\Xi^-$: reproduction of Fig. 4(a)")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def figure5(model: RPRModel, output: Path) -> None:
    energies = [1.95, 1.97, 2.07, 2.11, 2.14, 2.24, 2.28, 2.33, 2.42, 2.48, 2.79, 3.00]
    cth = np.linspace(-1.0, 1.0, 181)
    fig, axes = plt.subplots(4, 3, figsize=(10.0, 9.5), sharex=True, sharey=True)
    for ax, energy in zip(axes.flat, energies):
        u, sch, full = [], [], []
        for c in cth:
            kin = model.kinematics(energy, float(c))
            components = model.amplitude_components(energy, float(c))
            u.append(observables_from_components(model, kin, components, U_CHANNEL)[0])
            sch.append(observables_from_components(model, kin, components, S_CHANNEL)[0])
            full.append(observables_from_components(model, kin, components, None)[0])
        ax.plot(cth, u, "--", color="#d55e00", lw=1.1)
        ax.plot(cth, sch, ":", color="#0072b2", lw=1.2)
        ax.plot(cth, full, color="black", lw=1.5)
        ax.text(0.96, 0.88, f"{energy:.2f}", ha="right", transform=ax.transAxes)
        ax.grid(alpha=0.18)
        ax.set_ylim(0, 60)
    axes[0, 0].plot([], [], "--", color="#d55e00", label="u-channel Regge")
    axes[0, 0].plot([], [], ":", color="#0072b2", label="s-channel")
    axes[0, 0].plot([], [], color="black", label="full")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.supxlabel(r"$\cos\theta$")
    fig.supylabel(r"$d\sigma/d\Omega$ [$\mu$b/sr]")
    fig.suptitle(r"$K^-p\rightarrow K^+\Xi^-$: reproduction of Fig. 5(a)", y=0.995)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def figure7(model: RPRModel, output: Path) -> None:
    energies = [2.11, 2.24, 2.28, 2.42, 2.48]
    cth = np.linspace(-0.999, 0.999, 220)
    fig, axes = plt.subplots(len(energies), 1, figsize=(7.2, 9.0), sharex=True, sharey=True)
    for ax, energy in zip(axes, energies):
        u, sch, full = [], [], []
        for c in cth:
            kin = model.kinematics(energy, float(c))
            components = model.amplitude_components(energy, float(c))
            # The paper defines its plotted normal as x-hat cross z-hat,
            # opposite to the E42 convention K- cross K+ used by RPRModel.
            u.append(-observables_from_components(model, kin, components, U_CHANNEL)[1])
            sch.append(-observables_from_components(model, kin, components, S_CHANNEL)[1])
            full.append(-observables_from_components(model, kin, components, None)[1])
        ax.plot(cth, u, "--", color="#d55e00", lw=1.0)
        ax.plot(cth, sch, ":", color="#0072b2", lw=1.1)
        ax.plot(cth, full, color="black", lw=1.5)
        ax.axhline(0.0, color="0.65", lw=0.7)
        ax.text(0.96, 0.78, f"{energy:.2f} GeV", ha="right", transform=ax.transAxes)
        ax.grid(alpha=0.18)
        ax.set_ylim(-32, 25)
    axes[0].plot([], [], "--", color="#d55e00", label="u-channel Regge")
    axes[0].plot([], [], ":", color="#0072b2", label="s-channel")
    axes[0].plot([], [], color="black", label="full")
    axes[0].legend(frameon=False, fontsize=8, ncol=3)
    axes[-1].set_xlabel(r"$\cos\theta$")
    fig.supylabel(r"$P_y\,d\sigma/d\Omega$ [$\mu$b/sr]")
    fig.suptitle(r"$K^-p\rightarrow K^+\Xi^-$: reproduction of Fig. 7(a)", y=0.995)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    model = RPRModel()
    figure4(model, args.output_dir / "reproduce_fig4a.png")
    figure5(model, args.output_dir / "reproduce_fig5a.png")
    figure7(model, args.output_dir / "reproduce_fig7a.png")


if __name__ == "__main__":
    main()
