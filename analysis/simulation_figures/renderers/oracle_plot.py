"""Plot matched oracle and estimated nuisance assessments from validated data.

The RCT references are identical in the two nuisance-estimation settings.
"""

from pathlib import Path
import hashlib
import json
import string
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, FuncFormatter

ROOT = Path(__file__).resolve().parents[3]
DATA = Path(__file__).resolve().parents[1] / "inputs/oracle"
FIG = ROOT / "results/simulation-figures/figures"
QA = ROOT / "results/simulation-figures/validation"
STEM = "oracle_estimated_nuisance_comparison"
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
TITLES = [r"Cov. (θ)", r"Bias (θ)", r"MSE (θ)", "Width"]
GRID = np.arange(-12, 13) / 40


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    df = pd.read_csv(DATA / "summary.csv")
    protocol = json.loads((DATA / "protocol.json").read_text())
    assert protocol["validation"]["paired_datasets"] == 5000
    assert protocol["validation"]["max_direct_error"] < 1e-10
    assert len(df) == 200 and df.n_reps.eq(5000).all() and df.n_failed.eq(0).all()
    assert df.N_trial.eq(250).all() and np.isfinite(df[METRICS].to_numpy()).all()
    for method in ["ANOVA (RCT)", "ANCOVA (RCT)"]:
        top = df[(df["mode"] == MODES[0]) & (df.method == method)].sort_values("xi")
        bottom = df[(df["mode"] == MODES[1]) & (df.method == method)].sort_values("xi")
        assert np.array_equal(top[METRICS].to_numpy(), bottom[METRICS].to_numpy())

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )
    fig, axes = plt.subplots(2, 4, figsize=(6.75, 3.55), sharex=True, sharey="col")
    fig.subplots_adjust(
        left=0.105, right=0.987, bottom=0.22, top=0.94, wspace=0.42, hspace=0.57
    )
    curve_checks = []
    for i, mode in enumerate(MODES):
        for j, (metric, title) in enumerate(zip(METRICS, TITLES)):
            ax = axes[i, j]
            for method, label in METHODS:
                r = df[(df["mode"] == mode) & (df.method == method)].sort_values("xi")
                xx, yy = r.xi.to_numpy(), r[metric].to_numpy()
                assert len(r) == 25 and np.allclose(xx, GRID, atol=1e-14, rtol=0)
                if method.endswith("(RCT)"):
                    assert np.allclose(r.achieved_gap, 0, atol=1e-14, rtol=0)
                    assert np.ptp(yy) < 1e-12
                else:
                    assert np.allclose(r.achieved_gap, xx, atol=1e-12, rtol=0)
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
                fontsize=8.5,
                pad=4,
            )
            ax.set_xlim(-0.31, 0.31)
            ax.set_xticks([-0.3, 0, 0.3])
            ax.set_xticks([-0.15, 0.15], minor=True)

            # Label the endpoints and zero; retain the intermediate grid positions.
            def compact_tick(value, position):
                if np.isclose(value, 0):
                    return "0"
                label = f"{abs(value):.2f}"
                if abs(value) < 1:
                    label = label[1:]
                return ("\u2212" if value < 0 else "") + label

            ax.xaxis.set_major_formatter(FuncFormatter(compact_tick))
            ax.grid(which="both", color="#E6E6E6", lw=0.35)
            ax.set_axisbelow(True)
            ax.tick_params(length=2.2, width=0.45, pad=2)
            vals = df[metric].to_numpy()
            if metric == "coverage_theta":
                ax.set_ylim(0, 1.015)
                ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
                ax.axhline(0.95, color="#888888", lw=0.5, ls=(0, (3, 3)), zorder=0)
            elif metric == "bias_theta":
                ax.set_ylim(-0.34, 0.34)
                ax.set_yticks([-0.3, -0.15, 0, 0.15, 0.3])
                ax.axhline(0, color="#888888", lw=0.5, ls=(0, (3, 3)), zorder=0)
            else:
                span = vals.max() - vals.min()
                ax.set_ylim(
                    0 if metric == "mse_theta" else max(0, vals.min() - 0.09 * span),
                    vals.max() + max(0.005, 0.08 * span),
                )
                ax.yaxis.set_major_locator(MaxNLocator(4))
            if i == 1:
                ax.set_xlabel(
                    r"Causal-gap $\xi_{\lambda,\delta}$", fontsize=8.5, labelpad=3
                )
        pos = axes[i, 0].get_position()
        fig.text(
            0.014,
            (pos.y0 + pos.y1) / 2,
            mode,
            rotation=90,
            ha="left",
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )

    handles = [
        Line2D([], [], color=STYLES[m][0], ls=STYLES[m][1], lw=0.75, label=label)
        for m, label in METHODS
    ]
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.017),
        ncol=4,
        frameon=False,
        fontsize=8,
        handlelength=2.3,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist in list(axes.flat) + fig.texts + [legend]:
        b = artist.get_tightbbox(renderer)
        assert (
            b.x0 >= 0
            and b.y0 >= 0
            and b.x1 <= fig.bbox.width
            and b.y1 <= fig.bbox.height
        )
    for row in axes:
        for left_ax, right_ax in zip(row[:-1], row[1:]):
            assert (
                left_ax.get_tightbbox(renderer).x1 < right_ax.get_tightbbox(renderer).x0
            )
    lb = legend.get_window_extent(renderer)
    assert lb.y1 < min(ax.xaxis.label.get_window_extent(renderer).y0 for ax in axes[1])
    assert np.ptp([t.get_window_extent(renderer).y0 for t in legend.get_texts()]) < 0.01
    FIG.mkdir(exist_ok=True, parents=True)
    QA.mkdir(exist_ok=True, parents=True)
    outputs = []
    for ext in ["pdf"]:
        path = FIG / f"{STEM}.{ext}"
        kwargs = (
            {"metadata": {"CreationDate": None, "ModDate": None}}
            if ext == "pdf"
            else {}
        )
        fig.savefig(path, dpi=240, **kwargs)
        outputs.append(path)
    plt.close(fig)
    summary = {
        "success": True,
        "panels": 8,
        "curves": len(curve_checks),
        "legend_rows": 1,
        "input_sha256": sha(DATA / "summary.csv"),
        "curve_checks": curve_checks,
        "outputs": {p.name: sha(p) for p in outputs},
        "figure_inches": [6.75, 3.55],
        "font_points": {"ticks": 8, "legend": 8, "titles_axes_rows": 8.5},
        "row_order": MODES,
        "same_sample_counts": [167, 83, 500],
        "same_datasets": 5000,
    }
    (QA / "plot_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in summary.items() if k != "curve_checks"}, indent=2)
    )


if __name__ == "__main__":
    main()
