"""Exact operating characteristics indexed by a prespecified causal gap.

Source discrepancies are calibrated separately for each borrowing design.
Trial-only curves are explicitly zero-gap reference values, not biased RCTs.
"""

import json
import string
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from planning_paths import output_dir
from single_row_legends import fit_single_row

CURVES = output_dir("curves")
TABLES = output_dir("tables")
FIG = output_dir("figures")
VALIDATION = output_dir("validation")
XI = np.arange(-12, 13) / 40
Z = norm.ppf(0.975)
METRICS = ["coverage_theta", "bias_theta", "mse_theta", "width"]
HEADERS = [r"Cov. (θ)", r"Bias (θ)", r"MSE (θ)", "Width"]
COLORS = ["#333333", "#D55E00", "#0072B2", "#009E73", "#AA4499"]
LINES = [(0, (1, 2)), "--", "-", "-.", (0, (4, 1, 1, 1))]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.55,
    }
)


DISPLAY_LABELS = {
    "RCT: zero-gap reference": "RCT (100, 50): zero-gap reference",
}


def exact(label, nt, nc, ne, sigma=1, grid=XI, theta_alt=0.4):
    internal = nc / (nc + ne)
    variance = sigma * sigma * (1 / nt + 1 / (nc + ne))
    se = np.sqrt(variance)
    benchmark = ne == 0
    gap = np.zeros_like(grid) if benchmark else np.asarray(grid)
    delta = np.full_like(grid, np.nan) if benchmark else gap / (1 - internal)
    return pd.DataFrame(
        dict(
            method=label,
            xi_axis=grid,
            achieved_gap=gap,
            source_delta=delta,
            internal_share=internal,
            NT=nt,
            NC=nc,
            N0=ne,
            sigma=sigma,
            zero_gap_benchmark=benchmark,
            variance=variance,
            coverage_theta=norm.cdf(Z - gap / se) - norm.cdf(-Z - gap / se),
            bias_theta=gap,
            mse_theta=variance + gap**2,
            width=2 * Z * se,
            null_rejection=norm.cdf(gap / se - Z),
            power=norm.cdf((theta_alt + gap) / se - Z),
        )
    )


def decorate(ax, title):
    ax.set_title(title, loc="left", fontweight="bold", fontsize=7.1, pad=4)
    ax.grid(color="#E5E5E5", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=6.4, length=2.2, width=0.5, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#555555")
        spine.set_linewidth(0.5)
    ax.set_xlim(-0.31, 0.31)
    ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
    ax.set_xlabel(r"Causal gap $\xi_{\lambda,\delta}$", fontsize=7, labelpad=2)
    ax.set_box_aspect(1)


def save(fig, stem):
    fit_single_row(fig, stem)
    fig.savefig(FIG / f"{stem}.pdf")
    fig.savefig(FIG / f"{stem}.png", dpi=220)
    plt.close(fig)


def performance(rows, row_labels, stem):
    fig, axes = plt.subplots(2, 4, figsize=(6.75, 3.65))
    fig.subplots_adjust(
        left=0.073, right=0.988, bottom=0.18, top=0.955, wspace=0.38, hspace=0.48
    )
    methods = list(rows[0].method.unique())
    pooled = pd.concat(rows, ignore_index=True)
    for i, data in enumerate(rows):
        for j, (metric, name) in enumerate(zip(METRICS, HEADERS)):
            ax = axes[i, j]
            decorate(ax, f"({string.ascii_lowercase[i*4+j]}) {name}")
            if j == 0:
                ax.set_ylabel(row_labels[i], fontsize=6.5, labelpad=5)
            for k, method in enumerate(methods):
                d = data[data.method == method]
                ax.plot(
                    d.xi_axis,
                    d[metric],
                    color=COLORS[k],
                    ls=LINES[k],
                    lw=0.65,
                    marker="None",
                )
            if metric == "coverage_theta":
                ax.set_ylim(0, 1.02)
                ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
                ax.axhline(0.95, color="#666666", ls=(0, (2, 2)), lw=0.5, zorder=0)
            elif metric == "bias_theta":
                ax.set_ylim(-0.325, 0.325)
                ax.set_yticks([-0.3, 0, 0.3])
                ax.axhline(0, color="#666666", lw=0.4, zorder=0)
            else:
                ax.set_ylim(0, pooled[metric].max() * 1.07)
                ax.yaxis.set_major_locator(MaxNLocator(4))
    handles = [
        Line2D([], [], color=COLORS[k], ls=LINES[k], lw=0.65, marker="None", label=DISPLAY_LABELS.get(m, m))
        for k, m in enumerate(methods)
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.006),
        ncol=3 if len(methods) > 4 else 4,
        frameon=False,
        fontsize=6.3,
        handlelength=2.4,
        columnspacing=1.3,
    )
    save(fig, stem)


def main_figure():
    grid = np.linspace(-0.3, 0.3, 401)
    settings = [
        ("RCT: zero-gap reference", 100, 50, 0),
        ("ECT (100, 0)", 100, 0, 150),
        ("HCT (100, 50)", 100, 50, 150),
        ("HCT (200, 100)", 200, 100, 150),
        ("HCT (300, 150)", 300, 150, 150),
    ]
    data = pd.concat([exact(*s, grid=grid) for s in settings], ignore_index=True)
    data.to_csv(
        CURVES / "scalar_design_sensitivity.csv", index=False, float_format="%.16g"
    )
    fig, axes = plt.subplots(1, 3, figsize=(6.75, 2.8))
    fig.subplots_adjust(left=0.065, right=0.987, bottom=0.31, top=0.94, wspace=0.32)
    cols = ["coverage_theta", "null_rejection", "power"]
    titles = [
        r"(a) Cov. (θ)",
        "(b) Causal-null rejection",
        "(c) Efficacy declaration",
    ]
    for ax, col, title in zip(axes, cols, titles):
        decorate(ax, title)
        for k, (label, *_) in enumerate(settings):
            d = data[data.method == label]
            ax.plot(d.xi_axis, d[col], color=COLORS[k], ls=LINES[k], lw=0.75)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        if col == "coverage_theta":
            ax.axhline(0.95, color="#999999", ls=(0, (2, 2)), lw=0.45, zorder=0)
        elif col == "null_rejection":
            ax.axhline(0.025, color="#999999", ls=(0, (2, 2)), lw=0.45, zorder=0)
    handles = [
        Line2D([], [], color=COLORS[k], ls=LINES[k], lw=0.75, label=DISPLAY_LABELS.get(s[0], s[0]))
        for k, s in enumerate(settings)
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.015),
        ncol=3,
        frameon=False,
        fontsize=6.5,
        handlelength=2.5,
        columnspacing=1.5,
    )
    save(fig, "scalar_gap_design_sensitivity")
    assert np.allclose(data[data.zero_gap_benchmark].achieved_gap, 0)
    zero = data[np.isclose(data.xi_axis, 0)]
    assert np.allclose(zero.coverage_theta, 0.95) and np.allclose(
        zero.null_rejection, 0.025
    )
    return data


def reference_figures():
    signed = []
    scaled = []
    for ne in [150, 300]:
        settings = [
            ("RCT: zero-gap reference", 100, 50, 0),
            ("ECT (100, 0)", 100, 0, ne),
            ("HCT (100, 50)", 100, 50, ne),
            ("HCT (200, 100)", 200, 100, ne),
            ("HCT (300, 150)", 300, 150, ne),
        ]
        signed.append(
            pd.concat([exact(*s) for s in settings], ignore_index=True).assign(
                baseline_external_n=ne
            )
        )
        settings = [
            ("Baseline", 100, 50, ne, 1),
            ("Higher variance", 100, 50, ne, 1.25),
            ("All samples doubled", 200, 100, 2 * ne, 1),
            ("Trial doubled", 200, 100, ne, 1),
        ]
        scaled.append(
            pd.concat([exact(*s) for s in settings], ignore_index=True).assign(
                baseline_external_n=ne
            )
        )
    for stem, rows in [("reference_signed", signed), ("reference_scaling", scaled)]:
        df = pd.concat(rows, ignore_index=True)
        path = CURVES / f"scalar_{stem}.csv"
        df.to_csv(path, index=False, float_format="%.16g")
        plotted = pd.read_csv(path)
        rows = [plotted[plotted.baseline_external_n == n] for n in (150, 300)]
        prefix = "Baseline " if stem == "reference_scaling" else ""
        performance(
            rows,
            [prefix + r"$n_0=150$", prefix + r"$n_0=300$"],
            f"scalar_gap_{stem}_lines",
        )
        actual = df[~df.zero_gap_benchmark]
        assert np.allclose(
            (1 - actual.internal_share) * actual.source_delta, actual.xi_axis
        )
        assert np.allclose(actual.bias_theta, actual.xi_axis)
        assert np.allclose(actual.mse_theta, actual.variance + actual.xi_axis**2)


def precision_illustration():
    k = np.array([1, 4, 16, 64])
    se = np.sqrt(0.015 / k)
    half = Z * se
    fixed = np.full(4, 0.15)
    local = 0.15 / np.sqrt(k)
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.22, top=0.86, wspace=0.27)
    for j, (ax, gap, col) in enumerate(
        zip(axes, [fixed, local], ["#D55E00", "#7B3294"])
    ):
        title = "(a) Fixed causal gap" if j == 0 else "(b) Local causal gap"
        ax.set_title(title, loc="left", fontweight="bold", fontsize=8, pad=17)
        formula = (
            r"$\xi_{\lambda,\delta_k}=0.15$"
            if j == 0
            else r"$\xi_{\lambda,\delta_k}=0.15/\sqrt{k}$"
        )
        ax.text(
            0.5,
            1.018,
            formula,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=8,
        )
        ax.grid(color="#E5E5E5", lw=0.4)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color("#555555")
            spine.set_linewidth(0.55)
        y = np.arange(4)
        ax.errorbar(
            gap,
            y,
            xerr=half,
            fmt="o",
            color=col,
            lw=0.85,
            ms=2.4,
            capsize=3,
            capthick=0.65,
        )
        ax.plot(gap, y, color=col, ls="--", lw=0.6)
        for x, yy in zip(gap, y):
            ax.annotate(
                rf"$\xi_{{\lambda,\delta_k}}={x:g}$",
                (x, yy),
                xytext=(7, 9),
                textcoords="offset points",
                fontsize=6.8,
            )
        ax.axvline(0, color="#0072B2", ls=(0, (2, 2)), lw=0.8)
        ax.set_yticks(y, k)
        ax.invert_yaxis()
        ax.set_ylim(3.5, -0.5)
        ax.set_xlim(-0.15, 0.42)
        ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3, 0.4])
        ax.set_xlabel("Treatment-effect scale", fontsize=7.5)
        ax.set_ylabel(r"Sample-size multiplier $k$", fontsize=7.5)
        ax.tick_params(labelsize=7, length=2.5, width=0.5)
    handles = [
        Line2D(
            [],
            [],
            color="#0072B2",
            ls=(0, (2, 2)),
            lw=0.8,
            label=r"Trial effect: $\theta=0$",
        ),
        Line2D(
            [],
            [],
            color="#D55E00",
            marker="o",
            ms=2.5,
            lw=0.85,
            label="Fixed gap: TCE-centered span",
        ),
        Line2D(
            [],
            [],
            color="#7B3294",
            marker="o",
            ms=2.5,
            lw=0.85,
            label="Local gap: TCE-centered span",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.018),
        ncol=2,
        frameon=False,
        fontsize=7,
    )
    save(fig, "scalar_gap_precision_scaling")
    pd.DataFrame(
        dict(k=k, se=se, half_width=half, xi_fixed=fixed, xi_local=local)
    ).to_csv(TABLES / "precision_scaling.csv", index=False)
    assert np.allclose(local / se, 0.15 / np.sqrt(0.015))


if __name__ == "__main__":
    main_figure()
    reference_figures()
    precision_illustration()
    report = {
        "all_passed": True,
        "sensitivity_parameter": "population causal gap xi",
        "discrete_gap_grid": XI.tolist(),
        "source_calibration": "delta=xi/(1-lambda)",
        "RCT": "zero-gap reference values only; no nonzero RCT gap simulated",
        "analytical_not_Monte_Carlo": True,
    }
    (VALIDATION / "scalar_references.json").write_text(json.dumps(report, indent=2))
    print("Saved four scalar-gap reference figures.")
