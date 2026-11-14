"""Shared inputs and calculations for three fixed-Gaussian planning scenarios.

SEs are stipulated at their baseline values while borrowing gaps vary.
These reference calculations are not fitted-estimator simulation results.
"""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm


ALPHA = 0.025
THETA_STAR = 2.0
R_SQUARED = 0.25
GAMMA = 0.10
CRITICAL = float(norm.ppf(1 - ALPHA))
SCENARIOS = (
    {
        "id": "I",
        "sigma": 7.5,
        "n_treated": 120,
        "n_rct_control": 120,
        "n_hct_control": 60,
        "n_external": 300,
    },
    {
        "id": "II",
        "sigma": 10.0,
        "n_treated": 120,
        "n_rct_control": 120,
        "n_hct_control": 60,
        "n_external": 300,
    },
    {
        "id": "III",
        "sigma": 10.0,
        "n_treated": 240,
        "n_rct_control": 240,
        "n_hct_control": 120,
        "n_external": 600,
    },
)
STYLES = {
    "AIPW-U (HCT)": ("#176CB0", "-"),
    "AIPW (ECT)": ("#D56713", (0, (5, 2.5))),
    "ANOVA (RCT)": ("#525252", (0, (2, 2))),
    "ANCOVA (RCT)": ("#8250A5", (0, (6, 2, 1.5, 2))),
}
TEXT_SIZE = 10 * 72 / 72.27


def standard_errors(
    sigma, n_treated, n_rct_control, n_hct_control, n_external, r_squared
):
    v = sigma**2 * (1 - r_squared)
    return {
        "AIPW-U (HCT)": np.sqrt(v * (1 / n_treated + 1 / (n_hct_control + n_external))),
        "AIPW (ECT)": np.sqrt(v * (1 / n_treated + 1 / n_external)),
        "ANOVA (RCT)": sigma * np.sqrt(1 / n_treated + 1 / n_rct_control),
        "ANCOVA (RCT)": np.sqrt(v * (1 / n_treated + 1 / n_rct_control)),
    }


def scenario_ses(scenario):
    counts = tuple(
        scenario[k]
        for k in ("n_treated", "n_rct_control", "n_hct_control", "n_external")
    )
    assert np.isclose(counts[2] / (counts[2] + counts[3]), 1 / 6, rtol=0, atol=1e-15)
    return standard_errors(scenario["sigma"], *counts, R_SQUARED)


def probabilities(actual_gap, se, theta_star=THETA_STAR, critical=CRITICAL):
    """Coverage of theta, null rejection, and efficacy at theta_star."""
    standardized = np.asarray(actual_gap) / se
    return (
        norm.cdf(critical - standardized) - norm.cdf(-critical - standardized),
        norm.cdf(standardized - critical),
        norm.cdf(theta_star / se + standardized - critical),
    )


def is_borrowing(method):
    return "RCT" not in method


def threshold(se):
    return float(se * (CRITICAL + norm.ppf(GAMMA)))


def numeric_records():
    records = []
    for scenario in SCENARIOS:
        for method, se in scenario_ses(scenario).items():
            borrowing = is_borrowing(method)
            actual = np.array([-0.5, 0.0, 0.5]) if borrowing else np.zeros(3)
            coverage, rejection, efficacy = probabilities(actual, se)
            assert np.isclose(coverage[1], 0.95, atol=1e-14, rtol=0)
            assert np.isclose(rejection[1], ALPHA, atol=1e-14, rtol=0)
            assert np.isclose(
                efficacy[1], norm.cdf(THETA_STAR / se - CRITICAL), atol=1e-14
            )
            assert np.allclose(coverage, coverage[::-1], atol=1e-14, rtol=0)
            tol = threshold(se) if borrowing else None
            if borrowing:
                assert np.isclose(probabilities(tol, se)[1], GAMMA, atol=1e-14, rtol=0)
                assert probabilities(tol - 1e-6, se)[1] < GAMMA
                assert probabilities(tol + 1e-6, se)[1] > GAMMA
            else:
                assert all(
                    np.ptp(values) == 0 for values in (coverage, rejection, efficacy)
                )
            record = {
                **scenario,
                "method": method,
                "standard_error": float(se),
                "hct_weight": 1 / 6,
                "power_at_theta_star": float(efficacy[1]),
                "xi_tol_gamma_0_10": tol,
                "threshold_in_plotted_range": borrowing and tol <= 0.5,
                "curve_semantics": (
                    "matched_borrowing_gap" if borrowing else "zero_gap_reference"
                ),
            }
            for index, suffix in enumerate(("minus_0_5", "zero", "plus_0_5")):
                record[f"actual_gap_{suffix}"] = float(actual[index])
                record[f"coverage_{suffix}"] = float(coverage[index])
                record[f"rejection_{suffix}"] = float(rejection[index])
                record[f"efficacy_{suffix}"] = float(efficacy[index])
            records.append(record)
    return records


def write_numeric_report(root):
    tables = Path(root) / "tables"
    validation = Path(root) / "validation"
    tables.mkdir(parents=True, exist_ok=True)
    validation.mkdir(parents=True, exist_ok=True)
    records = numeric_records()
    report = {
        "calculation": "known_input_fixed_weight_Gaussian_reference",
        "not_fitted_estimator_simulation": True,
        "alpha_one_sided": ALPHA,
        "theta_star": THETA_STAR,
        "r_squared": R_SQUARED,
        "gamma": GAMMA,
        "gamma_note": "illustrative stress-test threshold",
        "hct_weight_all_scenarios": 1 / 6,
        "gap_calibration": "HCT: delta=6*xi/5; ECT: delta=xi; RCT: actual gap zero",
        "scenarios": list(SCENARIOS),
        "summary": records,
        "checks": {
            "coverage_at_zero": 0.95,
            "rejection_at_zero": ALPHA,
            "coverage_symmetric_in_gap": True,
            "threshold_equations_verified": True,
            "RCT_references_constant_at_actual_gap_zero": True,
        },
    }
    (validation / "planning_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    with (tables / "planning_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    return report


if __name__ == "__main__":
    from planning_paths import output_root

    write_numeric_report(output_root())
    print(output_root() / "tables" / "planning_summary.csv")
