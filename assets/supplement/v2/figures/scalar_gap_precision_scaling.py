"""Render the fixed-gap and local-gap Gaussian illustration."""
from pathlib import Path
import argparse
from matplotlib.font_manager import findfont
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "mathtext.fallback": "stixsans",
    "font.size": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.titleweight": "bold",
    "axes.linewidth": .55,
})

def precision_illustration(table, out):
    data = np.atleast_1d(np.genfromtxt(table, names=True, delimiter=","))
    k, half = data["k"], data["half_width"]
    fixed, local = data["xi_fixed"], data["xi_local"]
    if not np.array_equal(k, [1, 4, 16, 64]) or not np.all(np.isfinite([half, fixed, local])):
        raise ValueError("Expected four finite planning rows at k = 1, 4, 16, and 64.")
    findfont("Arial", fallback_to_default=False)
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.22, top=0.86, wspace=0.27)
    for j, (ax, gap, col) in enumerate(
        zip(axes, [fixed, local], ["#D55E00", "#7B3294"])
    ):
        title = "(a) Fixed causal gap" if j == 0 else "(b) Local causal gap"
        ax.set_title(title, loc="left", fontsize=8, pad=17)
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
        ax.set_yticks(y, [f"{value:g}" for value in k])
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
        bbox_to_anchor=(0.5, 0.018),
        ncol=3,
        frameon=False,
        fontsize=7,
    )
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "scalar_gap_precision_scaling.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path, help="precision_scaling.csv from run_planning.R")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    precision_illustration(args.table, args.output_dir)
