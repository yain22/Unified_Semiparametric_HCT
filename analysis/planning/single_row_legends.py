"""Fit analytical legends on one row without reducing their font sizes."""

import hashlib
import json

import numpy as np

from planning_paths import output_dir


def curve_state(fig):
    state = []
    for ax in fig.axes:
        lines = []
        for line in ax.lines:
            xy = np.column_stack([line.get_xdata(), line.get_ydata()])
            lines.append(hashlib.sha256(xy.tobytes()).hexdigest())
        state.append(
            dict(
                lines=lines,
                position=list(ax.get_position().bounds),
                xlim=list(ax.get_xlim()),
                ylim=list(ax.get_ylim()),
            )
        )
    return state


def fit_single_row(fig, stem):
    before = curve_state(fig)
    report = []
    for original in list(fig.legends):
        texts = original.get_texts()
        labels = [t.get_text() for t in texts]
        if original._ncols >= len(labels):
            report.append(dict(labels=labels, status="already_one_row"))
            continue
        fig.canvas.draw()
        anchor = fig.transFigure.inverted().transform(
            original.get_bbox_to_anchor().get_points()
        )[0]
        font = texts[0].get_fontproperties().copy()
        margin = 8 * fig.dpi / 72
        accepted = None
        candidates = [
            (original.columnspacing, original.handlelength),
            (1.0, 2.2),
            (0.8, 2.0),
            (0.6, 1.8),
        ]
        for spacing, length in candidates:
            candidate = fig.legend(
                handles=original.legend_handles,
                labels=labels,
                loc="lower center",
                bbox_to_anchor=(0.5, anchor[1]),
                ncol=len(labels),
                frameon=original.get_frame_on(),
                prop=font,
                handlelength=length,
                columnspacing=spacing,
                handletextpad=original.handletextpad,
                borderpad=original.borderpad,
                borderaxespad=original.borderaxespad,
            )
            fig.canvas.draw()
            bounds = candidate.get_window_extent(fig.canvas.get_renderer())
            if bounds.x0 >= margin and bounds.x1 <= fig.bbox.width - margin:
                accepted = dict(
                    labels=labels,
                    status="one_row",
                    fontsize=font.get_size_in_points(),
                    columnspacing=spacing,
                    handlelength=length,
                    horizontal_bounds_pixels=[bounds.x0, bounds.x1],
                    canvas_width_pixels=fig.bbox.width,
                )
                original.remove()
                break
            candidate.remove()
        report.append(
            accepted
            or dict(
                labels=labels,
                status="kept_for_legibility",
                fontsize=font.get_size_in_points(),
            )
        )
    assert before == curve_state(fig)
    validation = output_dir("validation") / "legends"
    validation.mkdir(parents=True, exist_ok=True)
    (validation / f"{stem}.json").write_text(
        json.dumps(dict(legends=report, curves_axes_unchanged=True), indent=2),
        encoding="utf-8",
    )
    return report
