"""Plot the four additional simulation scenarios."""

from pathlib import Path
import hashlib
import json
import string

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.special import roots_hermitenorm

ROOT = Path(__file__).resolve().parents[3]
DATA = Path(__file__).resolve().parents[1] / "inputs/extended"
FIG = ROOT / "results/simulation-figures/figures"
QA = ROOT / "results/simulation-figures/validation"
MODES = ["Oracle nuisance", "Estimated nuisance"]
METHODS = [
    ("U", "AIPW-U (HCT)"),
    ("AIPW (ECT)", "AIPW (ECT)"),
    ("ANOVA (RCT)", "ANOVA (RCT)"),
    ("ANCOVA (RCT)", "Linear ANCOVA (RCT)"),
]
STYLES = {
    "U": ("#2166AC", "-"),
    "AIPW (ECT)": ("#D95F02", "-."),
    "ANOVA (RCT)": ("#333333", (0, (1, 2))),
    "ANCOVA (RCT)": ("#7B4F9D", (0, (5, 2, 1, 2))),
}
METRICS = ["coverage_theta", "bias_theta", "mse_theta", "width"]
ALL_METRICS = METRICS + ["coverage_TCE", "bias_TCE", "estimated_variance"]
TITLES = [r"Cov. (θ)", r"Bias (θ)", r"MSE (θ)", "Width"]
GRID = np.arange(-12, 13) / 40


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    df = pd.read_csv(DATA / "summary.csv")
    protocol = json.loads((DATA / "protocol.json").read_text())
    validation = json.loads((DATA / "validation.json").read_text())
    calibration = pd.read_csv(DATA / "calibration.csv")
    assert validation["success"] and validation["failures"] == 0
    assert validation["paired_datasets"] == 20000
    assert len(df) == 800 and df.n_reps.eq(5000).all() and df.n_failed.eq(0).all()
    assert np.isfinite(df[ALL_METRICS].to_numpy()).all()
    gx, gw = roots_hermitenorm(128)
    gw /= np.sqrt(2 * np.pi)
    calibration_errors = {}
    for row in calibration.itertuples():
        odds_log = (
            np.log(row.nC / row.nE)
            + row.source_shift**2 / 2
            - (row.source_shift + 2 * row.variance_slope) * gx
        )
        external = 1 - expit(odds_log)
        shape = gx + row.quadratic_discrepancy * (gx**2 - 1)
        err = max(
            abs(np.sum(gw * external) - row.K0),
            abs(np.sum(gw * external * shape) - row.K1),
        )
        assert err < 1e-10
        calibration_errors[row.scenario] = float(err)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
        }
    )
    outputs = {}
    for scenario in df.scenario.drop_duplicates():
        selected = df[df.scenario.eq(scenario)].copy()
        raw_validation = {
            "mode": "rendering from previously validated summary; no simulation or reaggregation rerun"
        }
        fig, axes = plt.subplots(2, 4, figsize=(6.75, 3.55), sharex=True, sharey="col")
        fig.subplots_adjust(
            left=0.105, right=0.987, bottom=0.19, top=0.94, wspace=0.34, hspace=0.57
        )
        curve_checks = []
        for i, mode in enumerate(MODES):
            for j, (metric, title) in enumerate(zip(METRICS, TITLES)):
                ax = axes[i, j]
                for method, label in METHODS:
                    rows = selected[
                        selected["mode"].eq(mode) & selected.method.eq(method)
                    ].sort_values("xi")
                    xx, yy = rows.xi.to_numpy(), rows[metric].to_numpy()
                    assert len(rows) == 25 and np.allclose(xx, GRID, rtol=0, atol=1e-14)
                    if method.endswith("(RCT)"):
                        assert np.ptp(yy) < 1e-12 and rows.achieved_gap.eq(0).all()
                    color, ls = STYLES[method]
                    (line,) = ax.plot(xx, yy, color=color, ls=ls, lw=0.75)
                    assert np.array_equal(line.get_ydata(), yy)
                    curve_checks.append(
                        {
                            "mode": mode,
                            "method": method,
                            "metric": metric,
                            "xy_sha256": hashlib.sha256(
                                np.column_stack([xx, yy]).tobytes()
                            ).hexdigest(),
                        }
                    )
                ax.set_title(
                    f"({string.ascii_lowercase[4*i+j]}) {title}",
                    loc="left", fontweight="bold",
                    fontsize=7.1,
                    pad=4,
                )
                ax.set_xlim(-0.31, 0.31)
                ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
                ax.grid(color="#E6E6E6", lw=0.35)
                ax.set_axisbelow(True)
                ax.tick_params(length=2.2, width=0.45, pad=2)
                values = selected[metric].to_numpy()
                if metric == "coverage_theta":
                    ax.set_ylim(0, 1.015)
                    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
                    ax.axhline(0.95, color="#888888", lw=0.5, ls=(0, (3, 3)), zorder=0)
                elif metric == "bias_theta":
                    ax.set_ylim(-0.34, 0.34)
                    ax.set_yticks([-0.3, -0.15, 0, 0.15, 0.3])
                    ax.axhline(0, color="#888888", lw=0.5, ls=(0, (3, 3)), zorder=0)
                else:
                    span = values.max() - values.min()
                    ax.set_ylim(
                        (
                            0
                            if metric == "mse_theta"
                            else max(0, values.min() - 0.09 * span)
                        ),
                        values.max() + max(0.005, 0.08 * span),
                    )
                    ax.yaxis.set_major_locator(MaxNLocator(4))
                if i == 1:
                    ax.set_xlabel(
                        r"Causal-gap $\xi_{\lambda,\delta}$", fontsize=7, labelpad=3
                    )
            pos = axes[i, 0].get_position()
            fig.text(
                0.014,
                (pos.y0 + pos.y1) / 2,
                mode,
                rotation=90,
                ha="left",
                va="center",
                fontsize=7.2,
                fontweight="bold",
            )
        handles = [
            Line2D([], [], color=STYLES[m][0], ls=STYLES[m][1], lw=0.75, label=label)
            for m, label in METHODS
        ]
        legend = fig.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.54, 0.027),
            ncol=4,
            frameon=False,
            fontsize=6.2,
            handlelength=2.3,
            columnspacing=1.0,
            handletextpad=0.4,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for artist in list(axes.flat) + fig.texts + [legend]:
            bounds = artist.get_tightbbox(renderer)
            assert bounds.x0 >= 0 and bounds.y0 >= 0
            assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
        assert legend.get_window_extent(renderer).y1 < min(
            ax.xaxis.label.get_window_extent(renderer).y0 for ax in axes[1]
        )
        assert (
            np.ptp([t.get_window_extent(renderer).y0 for t in legend.get_texts()])
            < 0.01
        )
        for ax in axes.flat:
            for line in ax.lines:
                assert min(line.get_ydata()) >= ax.get_ylim()[0]
                assert max(line.get_ydata()) <= ax.get_ylim()[1]
        paths = []
        for ext in ["pdf"]:
            path = FIG / f"supp_scenario_{scenario}.{ext}"
            kwargs = (
                {"metadata": {"CreationDate": None, "ModDate": None}}
                if ext == "pdf"
                else {}
            )
            fig.savefig(path, dpi=300, **kwargs)
            paths.append(path)
        plt.close(fig)
        outputs[scenario] = {
            "raw_validation": raw_validation,
            "panels": 8,
            "legend_rows": 1,
            "curve_checks": curve_checks,
            "outputs": {p.name: sha(p) for p in paths},
        }
        print(f"Validated and plotted {scenario}", flush=True)
    report = {
        "success": True,
        "summary_sha256": sha(DATA / "summary.csv"),
        "calibration_independent_python_errors": calibration_errors,
        "validation_scope": "plot coordinates checked against original renderer hashes",
        "outputs": outputs,
    }
    (QA / "extended_plot_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
