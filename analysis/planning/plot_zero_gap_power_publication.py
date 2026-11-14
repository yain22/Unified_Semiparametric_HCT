"""SIM-style baseline power for the three shared study-design scenarios."""

import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from scipy.stats import norm

from study_design_scenarios import (
    SCENARIOS,
    STYLES,
    TEXT_SIZE,
    THETA_STAR,
    CRITICAL,
    scenario_ses,
    write_numeric_report,
)


DISPLAY_LABELS = {
    "AIPW-U (HCT)": "HCT reference",
    "AIPW (ECT)": "ECT reference",
}


def main():
    from planning_paths import output_root

    root = output_root()
    report = write_numeric_report(root)
    output = root / "figures"
    curve_dir = root / "curves"
    validation = root / "validation"
    for directory in (output, curve_dir, validation):
        directory.mkdir(parents=True, exist_ok=True)
    theta = np.linspace(-2, 4, 1201)
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
    fig, axes = plt.subplots(1, 3, figsize=(6.75, 2.95), sharey=True)
    fig.subplots_adjust(left=0.087, right=0.988, bottom=0.245, top=0.875, wspace=0.25)
    points, curves = [], []
    for panel, (ax, scenario) in enumerate(zip(axes, SCENARIOS)):
        for method, se in scenario_ses(scenario).items():
            color, style = STYLES[method]
            power_curve = norm.cdf(theta / se - CRITICAL)
            ax.plot(
                theta,
                power_curve,
                color=color,
                linestyle=style,
                linewidth=1.15,
                label=DISPLAY_LABELS.get(method, method),
            )
            power = float(norm.cdf(THETA_STAR / se - CRITICAL))
            open_marker = method == "AIPW (ECT)"
            ax.plot(
                THETA_STAR,
                power,
                marker="o",
                markersize=6.3 if open_marker else 3.5,
                markerfacecolor="none" if open_marker else color,
                markeredgecolor=color,
                markeredgewidth=0.85,
                linestyle="none",
                zorder=5,
            )
            expected = next(
                r
                for r in report["summary"]
                if r["id"] == scenario["id"] and r["method"] == method
            )
            assert np.isclose(
                power, expected["power_at_theta_star"], rtol=0, atol=1e-14
            )
            assert np.all(np.diff(power_curve) >= 0)
            assert np.all((power_curve >= 0) & (power_curve <= 1))
            points.append(
                {
                    "scenario": scenario["id"],
                    "sigma": scenario["sigma"],
                    "method": method,
                    "standard_error": float(se),
                    "theta": THETA_STAR,
                    "power": power,
                }
            )
            curves.extend(
                {
                    "scenario": scenario["id"],
                    "method": method,
                    "theta": float(x),
                    "power": float(y),
                }
                for x, y in zip(theta, power_curve)
            )
        ax.axhline(
            0.8, color="#999999", linewidth=0.65, linestyle=(0, (4, 3)), zorder=1
        )
        ax.axvline(
            THETA_STAR, color="black", linewidth=0.9, linestyle=(0, (1.5, 3)), zorder=2
        )
        ax.set(
            xlim=(-2.15, 4.15),
            ylim=(-0.04, 1.04),
            xlabel=r"$\theta$ (treatment effect)",
        )
        if panel == 0:
            ax.set_ylabel("Power", labelpad=4)
        ax.set_title(
            f"({chr(97+panel)}) σ = {scenario['sigma']:g}", loc="left", fontweight="bold", pad=5
        )
        ax.set_xticks([-2, 0, 2, 4])
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params(length=2.5, width=0.65, pad=2)
        ax.tick_params(which="minor", length=0)
        ax.xaxis.labelpad = 3
        ax.grid(which="major", color="#DFDFDF", linewidth=0.55)
        ax.grid(which="minor", color="#EEEEEE", linewidth=0.35)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color("#333333")
            spine.set_linewidth(0.65)
    handles = [
        Line2D(
            [],
            [],
            color=color,
            linestyle=style,
            linewidth=1.15,
            marker="o",
            markersize=5 if method == "AIPW (ECT)" else 3.5,
            markerfacecolor="none" if method == "AIPW (ECT)" else color,
            markeredgewidth=0.8,
            label=DISPLAY_LABELS.get(method, method),
        )
        for method, (color, style) in STYLES.items()
    ]
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.015),
        ncol=4,
        frameon=False,
        handlelength=2.1,
        columnspacing=0.95,
        handletextpad=0.45,
        borderaxespad=0,
        fontsize=TEXT_SIZE,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    checks = []
    for ax in axes:
        tight = ax.get_tightbbox(renderer)
        assert tight.x0 >= 0 and tight.y0 >= 0
        assert tight.x1 <= fig.bbox.width and tight.y1 <= fig.bbox.height
        title = ax._left_title.get_window_extent(renderer)
        bounds = ax.get_window_extent(renderer)
        assert title.x0 >= bounds.x0 - 0.5 and title.x1 <= bounds.x1 + 0.5
        checks.append(list(tight.bounds))
    lb = legend.get_window_extent(renderer)
    assert lb.x0 >= 0 and lb.x1 <= fig.bbox.width and lb.y0 >= 0
    assert lb.y1 < min(ax.get_tightbbox(renderer).y0 for ax in axes)
    for left, right in zip(axes[:-1], axes[1:]):
        assert left.get_tightbbox(renderer).x1 < right.get_tightbbox(renderer).x0
    for ext in ("pdf", "png", "svg"):
        fig.savefig(output / f"zero_gap_power.{ext}", dpi=360, facecolor="white")
    plt.close(fig)
    with (curve_dir / "baseline_power.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=curves[0].keys())
        writer.writeheader()
        writer.writerows(curves)
    (validation / "baseline_power.json").write_text(
        json.dumps(
            {
                "success": True,
                "scenarios": [r["id"] for r in SCENARIOS],
                "points": points,
                "figure_inches": [6.75, 2.95],
                "font_size_pdf_pt": TEXT_SIZE,
                "axes_bounds_pixels": checks,
                "legend_rows": 1,
                "artist_bounds_and_nonoverlap_verified": True,
            },
            indent=2,
        )
        + "\n"
    )
    print(output / "zero_gap_power.pdf")
    print("Three-scenario power curves and final-width typography verified.")


if __name__ == "__main__":
    main()
