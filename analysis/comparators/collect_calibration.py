"""Collect simulation calibration from retained summaries; write compact TeX tables.

Run from any directory. Inputs are the unchanged archived simulation summaries.
The resulting CSVs are portable aggregated evidence; no person-level records.
"""

from pathlib import Path
import argparse
import csv, hashlib, json, math
from statistics import NormalDist

here = Path(__file__).resolve().parent
root = here.parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir", type=Path, default=root / "results/calibration")
parser.add_argument("--simulation-dir", type=Path, help="Parent of completed baseline, oracle, and extended runs")
args = parser.parse_args()
output = args.output_dir.resolve()
if output == here:
    raise SystemExit(
        "Choose an output directory separate from the published reference tables."
    )
output.mkdir(parents=True, exist_ok=True)
files = {
    "baseline": here / "inputs/baseline-summary.csv",
    "transportability": here / "inputs/transportability-summary.csv",
    "extended": here / "inputs/extended-summary.csv",
    "oracle": here / "inputs/oracle-summary.csv",
}
if args.simulation_dir:
    input_root = args.simulation_dir.resolve()
    files = {
        "baseline": input_root / "baseline/summary.csv",
        "transportability": input_root / "baseline/transportability_reference_summary.csv",
        "extended": input_root / "extended/summary.csv",
        "oracle": input_root / "oracle/oracle_summary.csv",
    }
else:
    input_root = here
z = NormalDist().inv_cdf(0.975)
rows = []
for tag, p in files.items():
    for r in csv.DictReader(p.open(encoding="utf-8-sig")):
        if r["method"] not in ("U", "AIPW (ECT)"):
            continue
        n = int(r["N_trial"])
        row = {
            "setting": r.get("scenario", tag),
            "mode": r.get(
                "mode", "Estimated nuisance" if tag != "oracle" else "Oracle nuisance"
            ),
            "N_trial": n,
            "method": r["method"],
            "delta_identically_zero": int(tag == "transportability"),
            "gap": float(r.get("xi", 0)),
            "repetitions": int(r["n_reps"]),
            "bias_TCE": float(r["bias_TCE"]),
            "bias_mcse": float(r["bias_TCE_mcse"]),
            "empirical_sd": math.sqrt(float(r["empirical_variance"])),
            "mean_se": float(r["width"]) / (2 * z),
            "coverage_TCE": float(r["coverage_TCE"]),
            "coverage_mcse": float(r["coverage_TCE_mcse"]),
        }
        if "se" in r:
            assert abs(row["mean_se"] - float(r["se"])) < 1e-12
        rows.append(row)
out = output / "own_tce_calibration.csv"
with out.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=rows[0])
    w.writeheader()
    w.writerows(rows)


def lookup(setting, n, method, gap=0):
    x = [
        r
        for r in rows
        if r["setting"] == setting
        and r["N_trial"] == n
        and r["method"] == method
        and r["mode"] == "Estimated nuisance"
        and abs(r["gap"] - gap) < 1e-12
    ]
    assert len(x) == 1, (setting, n, method, gap, len(x))
    return x[0]


counts = {r["repetitions"] for r in rows}
if len(counts) != 1:
    raise ValueError("The calibration inputs must use the same number of repetitions.")
repetitions = counts.pop()
lines = [
    r"\begin{table}[!htbp]",
    r"\centering",
    r"\caption{Calibration for the statistical target, based on " + f"{repetitions:,}" + r" repetitions per row. Bias is measured against the statistical target of each procedure, SD is the empirical standard deviation, and SE is the mean reported standard error. The pointwise transportability check ($\delta\equiv0$) is separate from heterogeneous discrepancies with zero aggregate causal-gap. Additional scenarios are shown at zero aggregate gap.}\label{tab:own-tce-calibration}",
    r"\small",
    r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{3.5cm}lrrrrr@{}}",
    r"\toprule",
    r"Setting & Procedure & $n_{\mathrm{RCT}}$ & $\xi_{\lambda,\delta}$ & Bias & SD / SE & Coverage\\",
    r"\midrule",
]
for label, setting, n, method, gap in [
    (r"Pointwise $\delta\equiv0$", "transportability", 250, "U", 0),
    ("", "transportability", 500, "U", 0),
    ("", "transportability", 250, "AIPW (ECT)", 0),
    ("", "transportability", 500, "AIPW (ECT)", 0),
    ("Heterogeneous baseline", "baseline", 250, "U", -0.3),
    ("", "baseline", 250, "U", 0),
    ("", "baseline", 250, "U", 0.3),
    ("", "baseline", 500, "U", -0.3),
    ("", "baseline", 500, "U", 0),
    ("", "baseline", 500, "U", 0.3),
    ("", "baseline", 250, "AIPW (ECT)", 0),
    ("", "baseline", 500, "AIPW (ECT)", 0),
    ("Fewer external controls", "fewer_external", 250, "U", 0),
    ("Shifted source covariates", "source_shift", 250, "U", 0),
    ("Stronger heteroskedasticity", "stronger_heteroskedasticity", 250, "U", 0),
    ("Quadratic discrepancy", "quadratic_discrepancy", 250, "U", 0),
]:
    if label and len(lines) > 8:
        lines.append(r"\addlinespace[3pt]")
    r = lookup(setting, n, method, gap)
    name = "AIPW-U" if method == "U" else "ECT AIPW"
    lines.append(
        f"{label} & {name} & {n} & ${gap:.2f}$ & ${r['bias_TCE']:.4f}$ & ${r['empirical_sd']:.4f}$ / ${r['mean_se']:.4f}$ & ${r['coverage_TCE']:.4f}$ "
        + r"\\"
    )
lines.extend([r"\bottomrule", r"\end{tabular*}", r"\end{table}", ""])
(output / "own-tce-calibration-table.tex").write_text(
    "\n".join(lines), encoding="utf-8"
)
manifest = {
    "sources": {
        k: {
            "path": p.relative_to(input_root).as_posix(),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        for k, p in files.items()
    },
    "rows": len(rows),
    "mean_se": "Mean width/(2 z_.975), checked against saved mean SE where present",
    "distinction": "delta-identically-zero uses dedicated transportability_reference_summary; it is never inferred from zero aggregate gap.",
}
(output / "calibration-provenance.json").write_text(
    json.dumps(manifest, indent=2), encoding="utf-8"
)
print("Wrote calibration rows:", len(rows))
