"""Validate aggregate archives and plot the CGM estimates without participant data."""

from pathlib import Path
import argparse, json
from check_results import METHODS, validate_references
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
Z = 1.959963984540054
KEYS = METHODS


def guard_output(out, inputs):
    resolved = out.resolve()
    input_root = inputs.resolve()
    if resolved == input_root or input_root in resolved.parents:
        raise ValueError(
            "Output must be outside archived inputs, including their descendants."
        )
    return resolved


def validate(inputs):
    diagnostic = validate_references(inputs)
    rows = pd.read_csv(inputs / "primary/full_results.csv").iloc[:len(KEYS)]
    published = json.loads((inputs / "published_cgm_trials.json").read_text(encoding="utf-8-sig"))
    assert [r["study"] for r in published] == ["DIAMOND", "GOLD", "CITY"]
    return rows, published, diagnostic


def plot(rows, published, out):
    # Layout and geometry adapted from v3/code/focused_clinical_revision.py;
    # exact labels, values, axis positions, marker styles and dimensions retained.
    rows = rows.set_index("method").to_dict("index")
    names = [
        "ANOVA (RCT)",
        "ANCOVA (RCT)",
        "AIPW-R (HCT)",
        "AIPW-U (HCT)",
        "AIPW (ECT)",
    ]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(7.1, 4.2))
    fig.subplots_adjust(left=0.25, right=0.98, bottom=0.23, top=0.95)
    ticks = []
    labels = []
    values = []
    for i, (key, name) in enumerate(zip(KEYS, names)):
        yy = 7.6 - i
        row = rows[key]
        ticks.append(yy)
        labels.append(name)
        if i % 2 == 0:
            ax.axhspan(yy - 0.46, yy + 0.46, color="#F5F5F5", zorder=0)
        for kind, off, color, marker, style in [
            ("analytic", 0.12, "#0072B2", "o", "-"),
            ("bootstrap", -0.12, "#777777", "s", "--"),
        ]:
            lo, hi = row[kind + "_lower"], row[kind + "_upper"]
            line = ax.hlines(yy + off, lo, hi, color=color, ls=style, lw=1.4)
            point = ax.scatter(
                row["estimate"],
                yy + off,
                marker=marker,
                s=18,
                facecolor="white" if kind == "bootstrap" else color,
                edgecolor=color,
                zorder=3,
            )
            assert np.array_equal(
                line.get_segments()[0], [[lo, yy + off], [hi, yy + off]]
            )
            assert np.array_equal(point.get_offsets(), [[row["estimate"], yy + off]])
            values.append(
                dict(
                    method=key,
                    display_name=name,
                    interval=kind,
                    estimate=row["estimate"],
                    lower=lo,
                    upper=hi,
                )
            )
    heading = ax.text(
        -0.31,
        2.55,
        "Published randomized trials",
        transform=ax.get_yaxis_transform(),
        fontweight="bold",
        fontsize=8.5,
        ha="left",
        va="center",
        clip_on=False,
    )
    for i, row in enumerate(published):
        yy = 1.8 - i * 0.88
        ticks.append(yy)
        labels.append(row["study"] + " (" + row["author_year"].split("(")[1])
        lo, hi, estimate = row["ci_lower"], row["ci_upper"], row["estimate"]
        ax.axhspan(yy - 0.4, yy + 0.4, color="#F2F8F5", zorder=0)
        line = ax.hlines(yy, lo, hi, color="#009E73", lw=1.5)
        point = ax.scatter(estimate, yy, marker="D", color="#009E73", s=25, zorder=3)
        assert np.array_equal(line.get_segments()[0], [[lo, yy], [hi, yy]])
        assert np.array_equal(point.get_offsets(), [[estimate, yy]])
        values.append(
            dict(
                method=row["study"],
                interval="published",
                estimate=estimate,
                lower=lo,
                upper=hi,
            )
        )
    ax.set_yticks(ticks, labels)
    ax.set_xlim(-0.86, 0.22)
    ax.set_ylim(-0.5, 8.1)
    ax.set_xticks(np.arange(-0.8, 0.21, 0.2))
    ax.axvline(0, ls="--", lw=0.85, color="#777777")
    ax.grid(axis="x", color="#E6E6E6", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel(
        "HbA1c difference (percentage points; negative favors CGM)", labelpad=6
    )
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#888888")
    ax.tick_params(axis="y", length=0)
    handles = [
        Line2D(
            [],
            [],
            color="#0072B2",
            marker="o",
            markersize=4,
            lw=1.4,
            label="Working 95% interval",
        ),
        Line2D(
            [],
            [],
            color="#777777",
            marker="s",
            mfc="white",
            markersize=4,
            ls="--",
            lw=1.4,
            label="Bootstrap 95% interval",
        ),
        Line2D(
            [],
            [],
            color="#009E73",
            marker="D",
            markersize=4,
            lw=1.4,
            label="Published 95% interval",
        ),
    ]
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        frameon=False,
        ncol=3,
        fontsize=8.1,
        columnspacing=1.1,
        handletextpad=0.4,
        handlelength=1.7,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist in (
        ax.get_yticklabels() + ax.get_xticklabels() + [ax.xaxis.label, heading, legend]
    ):
        box = artist.get_window_extent(renderer)
        assert (
            box.x0 >= 0
            and box.x1 <= fig.bbox.width
            and box.y0 >= 0
            and box.y1 <= fig.bbox.height
        )
    assert (
        legend.get_window_extent(renderer).y1
        < ax.xaxis.label.get_window_extent(renderer).y0
    )
    assert len(values) == 13
    for suffix in ("pdf", "png"):
        fig.savefig(
            out / f"focused_clinical_comparison.{suffix}", dpi=220, facecolor="white"
        )
    plt.close(fig)
    (out / "figure_values.json").write_text(json.dumps(values, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "inputs")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    args.out = guard_output(args.out, args.inputs)
    args.out.mkdir(parents=True, exist_ok=True)
    rows, published, result = validate(args.inputs)
    if not args.validate_only:
        plot(rows, published, args.out)
    (args.out / "archive_validation.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
