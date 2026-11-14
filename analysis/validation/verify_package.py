"""Verify packaged code hashes and independently reaggregate saved numerical results.

Standard-library Python only. This does not refit patient data or rerun simulations.
For independent estimator checks, also run: Rscript analysis/paper_analysis.R validate
"""

from pathlib import Path
import argparse
import ast
import csv
import gzip
import hashlib
import json
import math
import statistics as st
import sys

ROOT = Path(__file__).resolve().parents[1]
NORMAL = st.NormalDist()
Z = NORMAL.inv_cdf(0.975)


def rows(relative):
    path = ROOT / relative
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT.parent / "results/package-check"
    )
    args = parser.parse_args()
    manifest = json.loads((ROOT / "code-manifest.json").read_text(encoding="utf-8"))
    for rel, expected in manifest["files"].items():
        path = ROOT / rel
        if (
            not path.is_file()
            or hashlib.sha256(path.read_bytes()).hexdigest() != expected
        ):
            raise ValueError("Missing or changed packaged file: " + rel)
        if path.suffix == ".py":
            ast.parse(path.read_text(encoding="utf-8-sig"), filename=rel)
    errors = []

    def same(actual, expected, name):
        error = abs(float(actual) - float(expected))
        if not math.isfinite(error) or error > 1e-11:
            raise ValueError(f"{name}: {actual} != {expected}")
        errors.append(error)

    records = rows("comparators/quadratic_rct_replicates.csv.gz")
    summaries = rows("comparators/quadratic_rct_summary.csv")
    for summary in summaries:
        group = [r for r in records if r["N_trial"] == summary["N_trial"]]
        if len(group) != int(summary["n_reps"]):
            raise ValueError("Quadratic repetition count mismatch")
        values = lambda field: [float(r[field]) for r in group]
        for target, field in [
            ("bias", "estimate"),
            ("mean_se", "se"),
            ("coverage", "coverage"),
            ("mse", "squared_error"),
            ("width", "width"),
        ]:
            same(st.mean(values(field)), summary[target], "quadratic " + target)
        same(
            st.stdev(values("estimate")),
            summary["empirical_sd"],
            "quadratic empirical SD",
        )
        for field in [
            "estimate",
            "coverage",
            "squared_error",
            "width",
            "paired_squared_error_difference",
            "paired_width_difference",
        ]:
            same(
                st.mean(values(field)),
                summary[field + "_mean"],
                "quadratic mean " + field,
            )
            same(
                st.stdev(values(field)) / math.sqrt(len(group)),
                summary[field + "_mcse"],
                "quadratic MCSE " + field,
            )

    sys.path.insert(0, str(ROOT / "clinical"))
    from check_results import validate_references
    clinical = validate_references(ROOT / "clinical/inputs")

    endpoints = rows("reference/common_discrepancy/endpoints.csv")
    for r in endpoints:
        variance = float(r["sampling_variance"])
        se = math.sqrt(variance)
        gap = float(r["source_discrepancy"]) * float(r["gap_multiplier"])
        same(se, r["standard_error"], "planning standard error")
        same(gap, r["causal_gap"], "planning gap")
        same(
            NORMAL.cdf(Z - gap / se) - NORMAL.cdf(-Z - gap / se),
            r["coverage"],
            "planning coverage",
        )
        same(NORMAL.cdf(gap / se - Z), r["null_rejection"], "planning null rejection")
        same(
            NORMAL.cdf((2 + gap) / se - Z),
            r["efficacy_declaration"],
            "planning efficacy",
        )
        same(variance + gap * gap, r["causal_mse"], "planning MSE")
        same(2 * Z * se, r["interval_width"], "planning interval width")

    report = {
        "passed": True,
        "hashed_files": len(manifest["files"]),
        "quadratic_replicates": len(records),
        "clinical": clinical,
        "planning_endpoint_rows": len(endpoints),
        "numerical_comparisons": len(errors),
        "maximum_absolute_error": max(errors),
        "scope": "Package integrity and independent saved-result reaggregation; no new fits.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "package-validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
