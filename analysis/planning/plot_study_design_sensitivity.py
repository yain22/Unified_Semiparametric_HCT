"""Signed-gap sensitivity for three shared fixed-Gaussian planning scenarios.

Fixed baseline SEs define an analytical planning reference, not fitted-estimator
simulation results. RCTs stay at actual gap zero throughout the plotted axis.
"""

import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np

from study_design_scenarios import (
    SCENARIOS,
    STYLES,
    TEXT_SIZE,
    ALPHA,
    GAMMA,
    THETA_STAR,
    CRITICAL,
    scenario_ses,
    is_borrowing,
    probabilities,
    threshold,
    write_numeric_report,
)


def main():
    from planning_paths import output_root

    root = output_root()
    report = write_numeric_report(root)
    output = root / "figures"
    curves = root / "curves"
    validation = root / "validation"
    for directory in (output, curves, validation):
        directory.mkdir(parents=True, exist_ok=True)
    gap = np.linspace(-0.5, 0.5, 2001)
    middle = len(gap) // 2
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": TEXT_SIZE,
            "axes.labelsize": TEXT_SIZE,
            "axes.titlesize": TEXT_SIZE,
            "xtick.labelsize": TEXT_SIZE,
            "ytick.labelsize": TEXT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(3, 3, figsize=(6.75, 6.8), sharex=True)
    fig.subplots_adjust(
        left=0.102, right=0.985, bottom=0.145, top=0.945, wspace=0.31, hspace=0.20
    )
    for ax, title in zip(
        axes[0], ("Causal coverage", "Null rejection", "Efficacy declaration")
    ):
        ax.set_title(title, loc="left", fontweight="bold", pad=8)
    rows, markers = [], []
    for row, scenario in enumerate(SCENARIOS):
        for method, se in scenario_ses(scenario).items():
            borrowing = is_borrowing(method)
            actual_gap = gap if borrowing else np.zeros_like(gap)
            values = probabilities(actual_gap, se)
            coverage, rejection, declaration = values
            assert np.isclose(coverage[middle], 0.95, rtol=0, atol=1e-14)
            assert np.isclose(rejection[middle], ALPHA, rtol=0, atol=1e-14)
            assert np.allclose(coverage, coverage[::-1], rtol=0, atol=1e-14)
            expected = next(
                r
                for r in report["summary"]
                if r["id"] == scenario["id"] and r["method"] == method
            )
            for index, suffix in ((0, "minus_0_5"), (middle, "zero"), (-1, "plus_0_5")):
                for key, curve in zip(("coverage", "rejection", "efficacy"), values):
                    assert np.isclose(
                        curve[index], expected[f"{key}_{suffix}"], rtol=0, atol=1e-14
                    )
            assert np.all(np.diff(rejection) >= -1e-14)
            assert np.all(np.diff(declaration) >= -1e-14)
            assert np.all(np.diff(coverage[: middle + 1]) >= -1e-14)
            assert np.all(np.diff(coverage[middle:]) <= 1e-14)
            if not borrowing:
                assert all(np.ptp(curve) == 0 for curve in values)
            color, style = STYLES[method]
            for col, curve in enumerate(values):
                assert np.all((curve >= 0) & (curve <= 1))
                (line,) = axes[row, col].plot(
                    gap,
                    curve,
                    color=color,
                    linestyle=style,
                    linewidth=1.1,
                    label=method,
                    zorder=3,
                )
                assert np.array_equal(line.get_ydata(), curve)
            tol = threshold(se) if borrowing else None
            if borrowing and gap[0] <= tol <= gap[-1]:
                ax = axes[row, 1]
                ax.vlines(
                    tol,
                    0,
                    GAMMA,
                    color=color,
                    linewidth=0.65,
                    linestyle=(0, (2, 2)),
                    alpha=0.65,
                    zorder=2,
                )
                ect = method == "AIPW (ECT)"
                ax.plot(
                    tol,
                    GAMMA,
                    marker="o",
                    markersize=5.4 if ect else 3.0,
                    markerfacecolor="none" if ect else color,
                    markeredgecolor=color,
                    markeredgewidth=0.8,
                    linestyle="none",
                    zorder=5,
                )
                markers.append(
                    {
                        "scenario": scenario["id"],
                        "method": method,
                        "gap": tol,
                        "null_rejection": GAMMA,
                    }
                )
            for x, actual, c, r, p in zip(gap, actual_gap, *values):
                rows.append(
                    {
                        "scenario": scenario["id"],
                        "method": method,
                        "plotted_gap": float(x),
                        "actual_gap": float(actual),
                        "standard_error": float(se),
                        "coverage": float(c),
                        "null_rejection": float(r),
                        "efficacy_declaration": float(p),
                    }
                )

        for col, ax in enumerate(axes[row]):
            ax.set_xlim(-0.52, 0.52)
            ax.set_xticks([-0.5, 0, 0.5], ["\u2212.5", "0", ".5"])
            ax.xaxis.set_minor_locator(MultipleLocator(0.25))
            ax.grid(which="major", color="#DFDFDF", linewidth=0.50)
            ax.grid(which="minor", axis="x", color="#EEEEEE", linewidth=0.35)
            ax.set_axisbelow(True)
            ax.tick_params(length=2.5, width=0.6, pad=2)
            ax.tick_params(which="minor", length=0)
            for spine in ax.spines.values():
                spine.set_color("#444444")
                spine.set_linewidth(0.6)
            if col == 0:
                ax.set_ylim(0.82, 0.965)
                ax.set_yticks([0.85, 0.90, 0.95], [".85", ".90", ".95"])
                guides = (0.95,)
            elif col == 1:
                ax.set_ylim(0, 0.20)
                ax.set_yticks(
                    [0, 0.05, 0.10, 0.15, 0.20], ["0", ".05", ".10", ".15", ".20"]
                )
                guides = (ALPHA, GAMMA)
            else:
                ax.set_ylim(0.30, 1.025)
                ax.set_yticks([0.4, 0.6, 0.8, 1.0], [".4", ".6", ".8", "1"])
                guides = (0.80,)
            for guide in guides:
                ax.axhline(
                    guide,
                    color="#999999",
                    linewidth=0.65,
                    linestyle=(0, (4, 3)),
                    zorder=1,
                )
        bounds = axes[row, 0].get_position()
        fig.text(
            0.018,
            (bounds.y0 + bounds.y1) / 2,
            scenario["id"],
            ha="center",
            va="center",
            fontsize=TEXT_SIZE,
        )

    xlabel = fig.text(
        0.545,
        0.083,
        r"Posited causal-gap $\xi_{\lambda,\delta}$",
        ha="center",
        fontsize=TEXT_SIZE,
    )
    labels = {"AIPW-U (HCT)": "HCT reference", "AIPW (ECT)": "ECT reference"}
    handles = [
        Line2D(
            [],
            [],
            color=color,
            linestyle=style,
            linewidth=1.1,
            label=labels.get(method, method),
        )
        for method, (color, style) in STYLES.items()
    ]
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.028),
        ncol=4,
        frameon=False,
        fontsize=TEXT_SIZE,
        handlelength=2.0,
        columnspacing=0.9,
        handletextpad=0.45,
        borderaxespad=0,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds_checks = []
    for ax in axes.flat:
        bounds = ax.get_tightbbox(renderer)
        assert bounds.x0 >= 0 and bounds.y0 >= 0
        assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
        bounds_checks.append(list(bounds.bounds))
        bottom, top = ax.get_ylim()
        for line in ax.lines:
            curve = np.asarray(line.get_ydata())
            assert np.all(curve >= bottom) and np.all(curve <= top)
        if ax.get_title():
            title = ax.title.get_window_extent(renderer)
            plot = ax.get_window_extent(renderer)
            assert title.x0 >= plot.x0 - 1 and title.x1 <= plot.x1 + 1
    for artists in axes:
        for left, right in zip(artists[:-1], artists[1:]):
            assert left.get_tightbbox(renderer).x1 < right.get_tightbbox(renderer).x0
    for row in range(2):
        assert min(ax.get_tightbbox(renderer).y0 for ax in axes[row]) > max(
            ax.get_tightbbox(renderer).y1 for ax in axes[row + 1]
        )
    for artist in [*fig.texts, legend]:
        bounds = artist.get_window_extent(renderer)
        assert bounds.x0 >= 0 and bounds.y0 >= 0
        assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
    assert legend.get_window_extent(renderer).y1 < xlabel.get_window_extent(renderer).y0
    assert xlabel.get_window_extent(renderer).y1 < min(
        ax.get_tightbbox(renderer).y0 for ax in axes[-1]
    )
    assert len(markers) == 4 and not any(r["scenario"] == "II" for r in markers)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(
            output / f"study_design_sensitivity.{ext}", dpi=300, facecolor="white"
        )
    plt.close(fig)
    with (curves / "matched_gap.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "success": True,
        "figure_inches": [6.75, 6.8],
        "font_size_pdf_pt": TEXT_SIZE,
        "rows": [r["id"] for r in SCENARIOS],
        "gap_range": [-0.5, 0.5],
        "points_per_curve": len(gap),
        "threshold_markers": markers,
        "legend_rows": 1,
        "bounds_pixels": bounds_checks,
        "layout_and_curve_bounds_verified": True,
        "all_baseline_and_endpoint_values_match_shared_numeric_report": True,
    }
    (validation / "matched_gap.json").write_text(json.dumps(result, indent=2) + "\n")
    print(output / "study_design_sensitivity.pdf")
    print("Three-scenario signed-gap curves, threshold markers and layout verified.")


if __name__ == "__main__":
    main()
