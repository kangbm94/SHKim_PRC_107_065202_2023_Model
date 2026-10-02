#!/usr/bin/env python3
"""Fit RPR amplitude scales to E42 polarization and optional cross sections.

The E42 prediction in each bin is folded over the accepted event sample:

    P_bin = sum_i w_i (d sigma_i/d Omega) P_i
            / sum_i w_i (d sigma_i/d Omega).

This is required for a quasifree data set whose events span sqrt(s) and angle.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Sequence

import numpy as np
from scipy.optimize import least_squares

from rpr_model import GEV2_TO_MICROBARN, PAULI_Y, RPRModel


def read_table(path: Path) -> list[dict[str, str]]:
    raw = [line for line in path.read_text().splitlines() if line.strip()]
    if not raw:
        raise ValueError(f"empty table: {path}")
    # The existing E42/HBC tables use a leading '#' on every header field.
    # Preserve that first line as the header and discard only later comments.
    header = raw[0].lstrip().lstrip("#")
    lines = [header] + [line for line in raw[1:] if not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"empty table: {path}")
    delimiter = "\t" if "\t" in lines[0] else ","
    return list(csv.DictReader(lines, delimiter=delimiter))


def get(row: dict[str, str], *names: str, default: str | None = None) -> str:
    lowered = {key.strip().lstrip("#").lower(): value for key, value in row.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    if default is not None:
        return default
    raise KeyError(f"none of {names} found in columns {tuple(row)}")


@dataclass
class Event:
    bin_id: str
    sqrt_s: float
    cos_theta: float
    weight: float


@dataclass
class PolarizationPoint:
    bin_id: str
    value: float
    error: float


@dataclass
class CrossSectionPoint:
    sqrt_s: float
    cos_theta: float
    value: float
    error: float


def load_events(path: Path) -> list[Event]:
    return [
        Event(
            get(row, "bin", "bin_id"),
            float(get(row, "sqrt_s", "sqrts")),
            float(get(row, "cos_theta", "costheta", "cth")),
            float(get(row, "weight", "event_weight", default="1")),
        )
        for row in read_table(path)
    ]


def load_polarization(path: Path) -> list[PolarizationPoint]:
    points = []
    for row in read_table(path):
        stat = float(get(row, "stat_error", "stat", "error", "err"))
        syst = float(get(row, "syst_error", "syst", default="0"))
        points.append(
            PolarizationPoint(
                get(row, "bin", "bin_id"),
                float(get(row, "polarization", "p_xi", "pxi", "value")),
                float(np.hypot(stat, syst)),
            )
        )
    return points


def load_cross_sections(path: Path) -> list[CrossSectionPoint]:
    return [
        CrossSectionPoint(
            float(get(row, "sqrt_s", "sqrts")),
            float(get(row, "cos_theta", "costheta", "cth")),
            float(get(row, "dsigma_domega", "cross_section", "xsection", "value")),
            float(get(row, "xsectionerr", "error", "err", "stat_error")),
        )
        for row in read_table(path)
    ]


class CachedPoint:
    def __init__(self, model: RPRModel, sqrt_s: float, cos_theta: float):
        self.components = model.amplitude_components(sqrt_s, cos_theta)
        kin = model.kinematics(sqrt_s, cos_theta)
        p = model.par
        self.prefactor = (
            p.m_p * p.m_xi / (16.0 * np.pi**2 * kin.s)
            * kin.p_final / kin.p_initial * GEV2_TO_MICROBARN
        )

    def observables(self, scales: Dict[str, float]) -> tuple[float, float]:
        amplitude = sum(
            (scales.get(name, 1.0) * value for name, value in self.components.items()),
            np.zeros((2, 2), complex),
        )
        density = amplitude @ amplitude.conj().T
        trace = float(np.trace(density).real)
        dsigma = self.prefactor * 0.5 * trace
        polarization = 0.0 if trace == 0.0 else float(np.trace(density @ PAULI_Y).real / trace)
        return dsigma, polarization


class E42Fit:
    def __init__(
        self,
        model: RPRModel,
        events: Sequence[Event],
        polarization: Sequence[PolarizationPoint],
        cross_sections: Sequence[CrossSectionPoint],
        parameter_names: Sequence[str],
        fit_normalization: bool,
    ):
        self.parameter_names = list(parameter_names)
        self.fit_normalization = fit_normalization
        self.polarization = list(polarization)
        self.cross_sections = list(cross_sections)
        self.event_cache = [(event, CachedPoint(model, event.sqrt_s, event.cos_theta)) for event in events]
        self.cross_cache = [CachedPoint(model, point.sqrt_s, point.cos_theta) for point in cross_sections]

        known_bins = {event.bin_id for event in events}
        missing = {point.bin_id for point in polarization} - known_bins
        if missing:
            raise ValueError(f"polarization bins without accepted MC events: {sorted(missing)}")
        if fit_normalization and not cross_sections:
            raise ValueError("overall normalization cannot be fitted without cross-section data")

    def unpack(self, values: np.ndarray) -> tuple[Dict[str, float], float]:
        scales = dict(zip(self.parameter_names, values[: len(self.parameter_names)]))
        normalization = float(values[-1]) if self.fit_normalization else 1.0
        return scales, normalization

    def predict_polarization(self, scales: Dict[str, float]) -> Dict[str, float]:
        numerator: Dict[str, float] = {}
        denominator: Dict[str, float] = {}
        for event, cached in self.event_cache:
            dsigma, polarization = cached.observables(scales)
            numerator[event.bin_id] = numerator.get(event.bin_id, 0.0) + event.weight * dsigma * polarization
            denominator[event.bin_id] = denominator.get(event.bin_id, 0.0) + event.weight * dsigma
        return {key: numerator[key] / denominator[key] for key in numerator}

    def residuals(self, values: np.ndarray) -> np.ndarray:
        scales, normalization = self.unpack(values)
        folded = self.predict_polarization(scales)
        residuals = [(folded[p.bin_id] - p.value) / p.error for p in self.polarization]
        for point, cached in zip(self.cross_sections, self.cross_cache):
            predicted, _ = cached.observables(scales)
            residuals.append((normalization * predicted - point.value) / point.error)
        return np.asarray(residuals)

    def fit(self):
        count = len(self.parameter_names) + int(self.fit_normalization)
        initial = np.ones(count)
        lower = np.full(count, -5.0)
        upper = np.full(count, 5.0)
        if self.fit_normalization:
            lower[-1], upper[-1] = 0.0, 10.0
        result = least_squares(self.residuals, initial, bounds=(lower, upper), jac="3-point")
        residual = self.residuals(result.x)
        chi2 = float(residual @ residual)
        ndof = len(residual) - len(result.x)
        covariance = None
        if result.jac.shape[0] >= result.jac.shape[1]:
            covariance = np.linalg.pinv(result.jac.T @ result.jac)
            if ndof > 0:
                covariance *= chi2 / ndof
        return result, chi2, ndof, covariance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--polarization", type=Path, required=True)
    parser.add_argument("--cross-section", type=Path)
    parser.add_argument(
        "--parameters",
        nargs="+",
        default=["Lambda2100", "Sigma2030"],
        help="amplitude scales to fit (default: dominant Lambda2100 and Sigma2030)",
    )
    parser.add_argument("--fit-normalization", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("fit_output"))
    args = parser.parse_args()

    model = RPRModel()
    events = load_events(args.events)
    polarization = load_polarization(args.polarization)
    cross_sections = load_cross_sections(args.cross_section) if args.cross_section else []
    fit = E42Fit(model, events, polarization, cross_sections, args.parameters, args.fit_normalization)
    result, chi2, ndof, covariance = fit.fit()
    scales, normalization = fit.unpack(result.x)
    errors = np.sqrt(np.diag(covariance)) if covariance is not None else np.full(len(result.x), np.nan)

    parameter_result = {
        name: {"value": float(value), "error": float(error)}
        for name, value, error in zip(
            args.parameters + (["normalization"] if args.fit_normalization else []), result.x, errors
        )
    }
    summary = {
        "success": bool(result.success),
        "message": result.message,
        "chi2": chi2,
        "ndof": ndof,
        "parameters": parameter_result,
        "covariance": None if covariance is None else covariance.tolist(),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "fit_result.json").write_text(json.dumps(summary, indent=2) + "\n")

    folded = fit.predict_polarization(scales)
    with (args.output_dir / "polarization_predictions.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["bin", "measured", "error", "predicted", "pull"])
        for point in polarization:
            predicted = folded[point.bin_id]
            writer.writerow([point.bin_id, point.value, point.error, predicted, (predicted - point.value) / point.error])

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
