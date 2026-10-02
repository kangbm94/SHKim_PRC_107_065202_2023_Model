#!/usr/bin/env python3
"""Diagnostic comparison to plot-extracted E42 bin centers.

This is intentionally *not* the final E42 fit.  It evaluates the free-proton
model at the centers of the published quasifree bins.  The proper calculation
uses fit_e42.py and the accepted event sample to fold the prediction over the
actual sqrt(s) and cos(theta) population in each bin.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

from fit_e42 import CachedPoint
from rpr_model import RPRModel


def read_points(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    return [
        {
            "sqrt_s": float(row["sqrts_center_GeV"]),
            "cos_theta": float(row["cosTheta_center"]),
            "value": float(row["PXi"]),
            "error": float(row["PXi_plot_err"]),
        }
        for row in rows
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data",
        type=Path,
        default=Path(__file__).parent / "data" / "e42_plot_extracted.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    args = parser.parse_args()

    model = RPRModel()
    points = read_points(args.data)
    cached = [CachedPoint(model, p["sqrt_s"], p["cos_theta"]) for p in points]

    def predictions(phi: float | None) -> np.ndarray:
        scales: dict[str, float] = {}
        if phi is not None:
            # A one-dimensional diagnostic deformation of the two dominant
            # resonances.  The fixed norm prevents polarization-only data from
            # driving both amplitudes to arbitrarily large values.
            scales = {
                "Lambda2100": float(np.sqrt(2.0) * np.cos(phi)),
                "Sigma2030": float(np.sqrt(2.0) * np.sin(phi)),
            }
        return np.asarray([item.observables(scales)[1] for item in cached])

    measured = np.asarray([p["value"] for p in points])
    errors = np.asarray([p["error"] for p in points])
    baseline = predictions(None)

    def chi2(phi: float) -> float:
        residual = (predictions(phi) - measured) / errors
        return float(residual @ residual)

    result = minimize_scalar(chi2, bounds=(-np.pi, np.pi), method="bounded")
    fitted = predictions(float(result.x))
    lambda_scale = float(np.sqrt(2.0) * np.cos(result.x))
    sigma_scale = float(np.sqrt(2.0) * np.sin(result.x))
    baseline_chi2 = float(np.sum(((baseline - measured) / errors) ** 2))
    fit_chi2 = chi2(float(result.x))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "warning": "Diagnostic bin-center comparison only; not acceptance folded.",
        "data_points": len(points),
        "baseline": {"chi2": baseline_chi2, "ndof": len(points)},
        "constrained_dominant_resonance_fit": {
            "constraint": "Lambda2100_scale^2 + Sigma2030_scale^2 = 2",
            "Lambda2100_scale": lambda_scale,
            "Sigma2030_scale": sigma_scale,
            "chi2": fit_chi2,
            "ndof": len(points) - 1,
        },
    }
    (args.output_dir / "e42_center_diagnostic.json").write_text(json.dumps(summary, indent=2) + "\n")

    with (args.output_dir / "e42_center_predictions.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["sqrt_s", "cos_theta", "measured", "error", "baseline", "constrained_fit"])
        for point, base, fit in zip(points, baseline, fitted):
            writer.writerow([point["sqrt_s"], point["cos_theta"], point["value"], point["error"], base, fit])

    angles = sorted({p["cos_theta"] for p in points})
    fig, axes = plt.subplots(1, len(angles), figsize=(11.0, 3.8), sharey=True)
    for ax, angle in zip(axes, angles):
        idx = [i for i, p in enumerate(points) if p["cos_theta"] == angle]
        energy = np.asarray([points[i]["sqrt_s"] for i in idx])
        order = np.argsort(energy)
        idx = np.asarray(idx)[order]
        energy = energy[order]
        ax.errorbar(energy, measured[idx], yerr=errors[idx], fmt="o", color="black", label="E42 plot extraction")
        ax.plot(energy, baseline[idx], "--", color="#d55e00", label="published parameters")
        ax.plot(energy, fitted[idx], "-", color="#0072b2", label="constrained diagnostic fit")
        ax.axhline(0.0, color="0.65", lw=0.7)
        ax.set_title(rf"$\cos\theta={angle:.3f}$")
        ax.set_xlabel(r"$\sqrt{s}$ [GeV]")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel(r"$P_y$ along $\mathbf{k}_{K^-}\times\mathbf{k}_{K^+}$")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("E42 bin-center diagnostic (free proton; no acceptance folding)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "e42_center_diagnostic.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
