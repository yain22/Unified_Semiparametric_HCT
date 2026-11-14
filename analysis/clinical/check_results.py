"""Check the published CGM summaries; optionally reaggregate local bootstrap refits."""

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist, stdev


METHODS = ["ANOVA (RCT)", "ANCOVA (RCT)", "AIPW-R", "AIPW-U (est)", "AIPW (ECT)"]
DIAGNOSTICS = ["matching_gap", "trial_augmented", "control_difference", "plugin_gap"]
Z = NormalDist().inv_cdf(0.975)


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_references(inputs, draws_root=None):
    errors = []
    counts = {}
    trial_points = None

    def close(actual, expected):
        difference = abs(float(actual) - float(expected))
        if not math.isfinite(difference) or difference > 1e-11:
            raise ValueError(f"Inconsistent CGM summary: {actual} versus {expected}")
        errors.append(difference)

    for scenario, n_external in (("primary", 273), ("subsequent", 256), ("near-index", 228)):
        data = inputs / scenario
        results = rows(data / "full_results.csv")
        assert [r["method"] for r in results] == METHODS + DIAGNOSTICS
        for r, point in zip(results, rows(data / "point_estimates.csv"), strict=True):
            name, estimate = r["method"], float(r["estimate"])
            assert point["method"] == name
            close(point["estimate"], estimate)
            se = float(r["bootstrap_se"])
            assert se > 0
            close(estimate - Z * se, r["bootstrap_lower"])
            close(estimate + Z * se, r["bootstrap_upper"])
            if name in METHODS:
                ase = float(r["analytic_se"])
                assert ase > 0
                close(ase, point["analytic_se"])
                close(estimate - Z * ase, r["analytic_lower"])
                close(estimate + Z * ase, r["analytic_upper"])
            else:
                assert all(r[f"analytic_{field}"] == "NA" for field in ("se", "lower", "upper"))
        estimates = {r["method"]: float(r["estimate"]) for r in results}
        close(estimates["AIPW-U (est)"], estimates["trial_augmented"] + estimates["matching_gap"])
        if trial_points is None:
            trial_points = [estimates[k] for k in METHODS[:2]]
        else:
            for k, initial in zip(METHODS[:2], trial_points):
                close(estimates[k], initial)
        cells = rows(data / "partition.csv")
        assert len(cells) == 15
        for group, expected in (("T", 229), ("C", 213), ("E", n_external)):
            close(sum(int(r["Freq"]) for r in cells if r["group"] == group), expected)
        expected_sizes = {"trial_treated": 228, "internal": 212, "external": n_external,
                          "ect_treated": 228, "ect_external": n_external}
        probabilities = rows(data / "source_probabilities.csv")
        assert {r["quantity"] for r in probabilities} == expected_sizes.keys()
        for r in probabilities:
            assert int(r["n"]) == expected_sizes[r["quantity"]]
            q = [float(r[v]) for v in ("min", "p05", "median", "p95", "max")]
            assert sorted(q) == q and 0 < q[0] <= q[-1] < 1
        if draws_root is not None:
            draws = rows(draws_root / scenario / "bootstrap_draws.csv")
            assert [int(r["replicate"]) for r in draws] == list(range(1, 2001))
            for draw in draws:
                close(draw["AIPW-U (est)"], float(draw["trial_augmented"]) + float(draw["matching_gap"]))
            for r in results:
                close(stdev(float(d[r["method"]]) for d in draws), r["bootstrap_se"])
            counts[scenario] = len(draws)
    return {"estimates_per_cohort": len(METHODS + DIAGNOSTICS), "cohorts": 3,
            "local_bootstrap_draws_checked": counts, "comparisons": len(errors),
            "max_difference": max(errors)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path(__file__).resolve().parent / "inputs")
    parser.add_argument("--local-refits", type=Path, help="Directory containing the three 2,000-refit CGM runs")
    args = parser.parse_args()
    print(json.dumps(validate_references(args.inputs, args.local_refits), indent=2))
