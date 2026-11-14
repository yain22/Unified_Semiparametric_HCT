"""Draw the observed-data model regions and the changes of observation law."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path as PlotPath

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path.cwd())
    args = parser.parse_args()
    OUT = args.output_dir.resolve()
    plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "cm",
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "savefig.facecolor": "white"})
    from matplotlib.font_manager import findfont
    findfont("Arial", fallback_to_default=False)
    fig = plt.figure(figsize=(3.25, 2.12))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 3.25), ylim=(-.16, 2.45))
    ax.axis("off")
    ink, blue, green, orange = "#172331", "#2667A8", "#167565", "#AE5C17"


    def txt(x, y, content, size=8, color=ink, weight="normal", **kwargs):
        return ax.text(x, y, content, ha="center", va="center",
                       fontsize=max(6., .82 * size), color=color,
                       fontweight=weight, **kwargs)


    ax.add_patch(FancyBboxPatch((.035, .035), 3.18, 2.38,
        boxstyle="round,pad=0,rounding_size=.06", facecolor="#FCFCFD",
        edgecolor="#586575", linewidth=.85))
    txt(1.625, 2.31, "Unified causal-gap sensitivity framework", 8.5, weight="bold")
    ax.add_patch(Ellipse((1.625, 1.535), 3.02, 1.33,
                        facecolor="#EDF4FC", edgecolor=blue, lw=1.15))
    txt(1.625, 2., "Unrestricted HCT\nobserved-data model", 8., blue,
        weight="bold", linespacing=1.10)
    txt(1.625, 1.765, r"$\mathcal{M}_{\mathrm{HCT},U}$", 11, blue)
    ax.add_patch(Ellipse((1.625, 1.29), 2.40, .65,
                        facecolor="#DFF1EB", edgecolor=green, lw=1.0))
    txt(1.625, 1.455, "HCT under transportability", 7.4, green, weight="bold")
    txt(1.625, 1.28, r"$\mathcal{M}_{\mathrm{HCT,trans}}\subset\mathcal{M}_{\mathrm{HCT},U}$", 9.4, green)
    txt(1.625, 1.125, r"$m_{1,0}(x)=m_0(x)\ \Longrightarrow\ \xi_{\lambda,\delta}=0$", 8.1, green)
    for vertices, color in [([(1.04, .922), (.815, .922), (.815, .635)], "#566373"),
                            ([(2.21, .922), (2.435, .922), (2.435, .635)], orange)]:
        path = PlotPath(vertices, [PlotPath.MOVETO, PlotPath.LINETO, PlotPath.LINETO])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=8.5,
                                    lw=.9, color=color, joinstyle="miter"))
    txt(.36, .77, "Remove external\ncontrol records", 6.8, "#566373", linespacing=1.15).set_fontsize(5.3)
    txt(2.89, .77, "Remove internal\ncontrol records", 6.8, orange, linespacing=1.15).set_fontsize(5.3)
    txt(1.12, .755, r"$\lambda\equiv1$", 9.5, "#465363")
    txt(2.13, .755, r"$\lambda\equiv0$", 9.5, orange)
    for x, fc, ec in [(.10, "#F1F3F6", "#697789"), (1.72, "#FFF3E8", orange)]:
        ax.add_patch(FancyBboxPatch((x, .12), 1.43, .50,
            boxstyle="round,pad=0,rounding_size=.045", facecolor=fc, edgecolor=ec, lw=.85))
    txt(.815, .48, "RCT", 9, "#465363", weight="bold")
    txt(.815, .27, r"$\mathcal{M}_{\mathrm{RCT}}$", 11, "#465363")
    txt(2.435, .48, "ECT", 9, orange, weight="bold")
    txt(2.435, .27, r"$\mathcal{M}_{\mathrm{ECT}}$", 11, orange)
    txt(1.625, -.075, "Design reductions redefine the observation law", 7.3)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist in ax.texts:
        box = artist.get_window_extent(renderer)
        assert 0 <= box.x0 < box.x1 <= fig.bbox.width, artist.get_text()
        assert 0 <= box.y0 < box.y1 <= fig.bbox.height, artist.get_text()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "observed_model_relationships.pdf", pad_inches=0)
    fig.savefig(OUT / "observed_model_relationships.png", dpi=360, pad_inches=0)
    plt.close(fig)
    print(OUT / "observed_model_relationships.pdf")


if __name__ == '__main__':
    main()
