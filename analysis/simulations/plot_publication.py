"""Render the six publication simulation figures to an explicit output directory.

Uses the checked-in simulation_figures renderer and input files by relative path.
Only complete 5000-repetition summaries are accepted by publication renderers.
"""

from pathlib import Path
import argparse, importlib.util, json, sys
import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
LABELS = HERE.parent / "simulation_figures"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    for name in ["baseline", "oracle", "extended"]:
        p.add_argument(
            "--" + name,
            type=Path,
            help="Directory of complete regenerated results; defaults to retained publication summaries",
        )
    args = p.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    reports = {}
    for name, flag in [
        ("heterogeneous", "baseline"),
        ("oracle", "oracle"),
        ("extended", "extended"),
    ]:
        spec = importlib.util.spec_from_file_location(
            "publication_" + name, LABELS / "renderers" / f"{name}_plot.py"
        )
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        m.DATA = (getattr(args, flag) or LABELS / "inputs" / name).resolve()
        m.ROOT = out
        m.FIG = out
        m.QA = out / "validation"
        m.QA.mkdir(exist_ok=True)
        attr = "MAIN" if name == "heterogeneous" else "METHODS"
        setattr(
            m,
            attr,
            [
                (key, "Linear ANCOVA (RCT)" if key == "ANCOVA (RCT)" else label)
                for key, label in getattr(m, attr)
            ],
        )
        if name == "heterogeneous":
            m.STYLES = {
                key: (color, ls, "None", filled)
                for key, (color, ls, marker, filled) in m.STYLES.items()
            }
            reports[name] = m.draw(
                pd.read_csv(m.DATA / "summary.csv"),
                m.MAIN,
                "heterogeneous_gap_design_comparison_lines",
            )
        else:
            m.main()
            report = (
                "plot_validation.json"
                if name == "oracle"
                else "extended_plot_validation.json"
            )
            reports[name] = json.loads((m.QA / report).read_text())
    (out / "validation/plot_manifest.json").write_text(
        json.dumps(reports, indent=2), encoding="utf-8"
    )
    print("Rendered six figures to", out)


if __name__ == "__main__":
    main()
