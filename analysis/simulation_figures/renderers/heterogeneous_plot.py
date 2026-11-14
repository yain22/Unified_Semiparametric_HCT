"""Plot operating characteristics for the baseline simulation."""

from pathlib import Path
import hashlib, json, string
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[3]
DATA = Path(__file__).resolve().parents[1] / "inputs/heterogeneous"
FIG = ROOT / "results/simulation-figures/figures"
QA = ROOT / "results/simulation-figures/validation"
GRID = np.arange(-12, 13) / 40
METRICS = ["coverage_theta", "bias_theta", "mse_theta", "width"]
TITLES = [r"Cov. (θ)", r"Bias (θ)", r"MSE (θ)", "Width"]
MAIN = [
    ("U", "AIPW-U (HCT)"),
    ("AIPW (ECT)", "AIPW (ECT)"),
    ("ANOVA (RCT)", "ANOVA (RCT)"),
    ("ANCOVA (RCT)", "Linear ANCOVA (RCT)"),
]
PAIRED = [
    ("U", "AIPW-U"),
    ("R-matched", "AIPW-R (matched share)"),
    ("R-unit-same-law", "AIPW-R (unit ratio)"),
]
STYLES = {
    "U": ("#2166AC", "-", "None", True),
    "AIPW (ECT)": ("#D95F02", "-.", "None", False),
    "ANOVA (RCT)": ("#333333", (0, (1, 2)), "None", True),
    "ANCOVA (RCT)": ("#7B4F9D", (0, (5, 2, 1, 2)), "None", True),
    "R-matched": ("#00876C", "--", "None", False),
    "R-unit-same-law": ("#2166AC", "--", "None", False),
}
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "axes.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def draw(data, methods, stem, same_law=False):
    selected = data[data.method.isin([m for m, _ in methods])].copy()
    assert selected.n_reps.eq(5000).all()
    if "n_failed" in selected:
        assert selected.n_failed.eq(0).all()
    assert np.isfinite(selected[METRICS].to_numpy()).all()
    fig, axes = plt.subplots(2, 4, figsize=(6.75, 3.7), sharex=True, sharey="col")
    fig.subplots_adjust(
        left=0.074, right=0.984, bottom=0.215, top=0.943, wspace=0.35, hspace=0.47
    )
    checks = []
    for i, n in enumerate([250, 500]):
        for j, (metric, title) in enumerate(zip(METRICS, TITLES)):
            ax = axes[i, j]
            for k, (method, label) in enumerate(methods):
                r = selected[
                    (selected.N_trial == n) & (selected.method == method)
                ].sort_values("xi_index")
                assert len(r) == 25
                xx = r.xi.to_numpy()
                yy = r[metric].to_numpy()
                assert np.allclose(xx, GRID, rtol=0, atol=1e-14)
                if same_law:
                    assert np.allclose(r.reference_gap, xx, rtol=0, atol=1e-14)
                elif method.endswith("(RCT)"):
                    assert np.allclose(r.achieved_gap, 0, rtol=0, atol=1e-14)
                else:
                    assert np.allclose(r.achieved_gap, xx, rtol=0, atol=1e-12)
                c, ls, mk, filled = STYLES[method]
                (line,) = ax.plot(
                    xx,
                    yy,
                    c=c,
                    ls=ls,
                    lw=0.67,
                    marker=mk,
                    ms=2.2,
                    mew=0.5,
                    markevery=(k % 5, 5),
                    mfc=c if filled else "white",
                )
                assert np.array_equal(line.get_xdata(), xx) and np.array_equal(
                    line.get_ydata(), yy
                )
                checks.append(
                    dict(
                        n=n,
                        method=method,
                        metric=metric,
                        xy_sha256=hashlib.sha256(
                            np.column_stack([xx, yy]).tobytes()
                        ).hexdigest(),
                    )
                )
            ax.set_title(
                f"({string.ascii_lowercase[i*4+j]}) {title}",
                loc="left", fontweight="bold",
                fontsize=7.1,
                pad=4,
            )
            ax.set_box_aspect(1)
            ax.set_xlim(-0.31, 0.31)
            ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
            if i == 1 and not same_law:
                ax.set_xticklabels(["\u2212.30", "", "0", "", ".30"])
            ax.grid(color="#E6E6E6", lw=0.35)
            ax.set_axisbelow(True)
            ax.tick_params(length=2.2, width=0.45, pad=2)
            vals = selected[metric].to_numpy()
            if metric == "coverage_theta":
                ax.set_ylim(0, 1.015)
                ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
                ax.axhline(0.95, c="#888888", lw=0.5, ls=(0, (3, 3)), zorder=0)
            elif metric == "bias_theta":
                bound = max(0.34, 1.05 * max(abs(vals)))
                ax.set_ylim(-bound, bound)
                ax.yaxis.set_major_locator(MaxNLocator(5))
                ax.axhline(0, c="#888888", lw=0.5, ls=(0, (3, 3)), zorder=0)
            else:
                span = vals.max() - vals.min()
                ax.set_ylim(
                    0 if metric == "mse_theta" else max(0, vals.min() - 0.08 * span),
                    vals.max() + max(0.005, 0.08 * span),
                )
                ax.yaxis.set_major_locator(MaxNLocator(4))
            if i == 1:
                ax.set_xlabel(
                    (
                        "Precision-share HCT gap\n" + r"$\xi_{\lambda,\delta}$"
                        if same_law
                        else r"Causal gap $\xi_{\lambda,\delta}$"
                    ),
                    fontsize=7,
                    labelpad=3,
                )
            if j == 0:
                ax.set_ylabel(rf"$n_{{\mathrm{{RCT}}}}={n}$", fontsize=7, labelpad=5)
    handles = []
    for m, label in methods:
        c, ls, mk, filled = STYLES[m]
        handles.append(
            Line2D(
                [],
                [],
                c=c,
                ls=ls,
                lw=0.67,
                marker=mk,
                ms=2.3,
                mew=0.5,
                mfc=c if filled else "white",
                label=label,
            )
        )
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.028),
        ncol=len(methods),
        frameon=False,
        fontsize=6.1,
        handlelength=2.1,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lb = legend.get_window_extent(renderer)
    assert lb.x0 >= 0 and lb.x1 <= fig.bbox.width and lb.y0 >= 0
    assert np.ptp([t.get_window_extent(renderer).y0 for t in legend.get_texts()]) < 0.01
    assert lb.y1 < min(ax.xaxis.label.get_window_extent(renderer).y0 for ax in axes[1])
    assert [t.get_text() for t in legend.get_texts()] == [label for _, label in methods]
    for ax in axes.flat:
        b = ax.get_window_extent(renderer)
        assert abs(b.width - b.height) < 1e-8
        b = ax.get_tightbbox(renderer)
        assert (
            b.x0 >= 0
            and b.x1 <= fig.bbox.width
            and b.y0 >= 0
            and b.y1 <= fig.bbox.height
        )
        for line in ax.lines:
            yy = line.get_ydata()
            assert min(yy) >= ax.get_ylim()[0] and max(yy) <= ax.get_ylim()[1]
    FIG.mkdir(exist_ok=True, parents=True)
    QA.mkdir(parents=True, exist_ok=True)
    outputs = []
    for ext in ["pdf"]:
        path = FIG / f"{stem}.{ext}"
        fig.savefig(
            path,
            dpi=240,
            **(
                {"metadata": {"CreationDate": None, "ModDate": None}}
                if ext == "pdf"
                else {}
            ),
        )
        outputs.append(path)
    plt.close(fig)
    return dict(
        methods=[label for _, label in methods],
        panels=8,
        curves=len(checks),
        legend_rows=1,
        axis=(
            "Precision-share HCT gap xi_lambda_delta"
            if same_law
            else "causal gap xi_lambda_delta"
        ),
        curve_checks=checks,
        outputs={p.name: sha(p) for p in outputs},
    )


def main():
    data = pd.read_csv(DATA / "summary.csv")
    result = draw(data, MAIN, "heterogeneous_gap_design_comparison_lines")
    (QA / "plot_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    main()
