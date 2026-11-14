"""Plot simulation results from the supplied aggregate summaries.

Run from any directory. Python requires numpy, pandas, matplotlib, and scipy.
Optional PDF geometry checks use a Python interpreter with PyMuPDF.
"""

from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys
import pandas as pd
from renderers import oracle_plot, heterogeneous_plot, extended_plot
import subprocess

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = HERE / "inputs"
REFERENCES = json.loads((INPUTS / "reference_curve_checks.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


RENDERERS = {
    "oracle": oracle_plot,
    "heterogeneous": heterogeneous_plot,
    "extended": extended_plot,
}


def render(destination, audit):
    destination.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)
    results = {}
    for name in ["oracle", "heterogeneous", "extended"]:
        module = RENDERERS[name]
        module.ROOT = ROOT
        module.DATA = INPUTS / name
        module.FIG = destination
        module.QA = audit
        if name == "heterogeneous":
            result = module.draw(
                pd.read_csv(module.DATA / "summary.csv"),
                module.MAIN,
                "heterogeneous_gap_design_comparison_lines",
            )
            assert result["curve_checks"] == REFERENCES[name]
            results[name] = result
        else:
            module.main()
            filename = (
                "plot_validation.json"
                if name == "oracle"
                else "extended_plot_validation.json"
            )
            result = json.loads((audit / filename).read_text())
            if name == "oracle":
                assert result["curve_checks"] == REFERENCES[name]
            else:
                for scenario, entry in result["outputs"].items():
                    assert entry["curve_checks"] == REFERENCES[name][scenario]
            results[name] = result
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "results/simulation-figures"
    )
    parser.add_argument(
        "--check-vectors",
        action="store_true",
        help="Also compare PDF geometry with the published assets (requires PyMuPDF)",
    )
    parser.add_argument(
        "--pdf-python",
        default=sys.executable,
        help="Python interpreter with PyMuPDF for --check-vectors",
    )
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs_before = {
        p.relative_to(INPUTS).as_posix(): sha(p)
        for p in INPUTS.rglob("*")
        if p.is_file()
    }
    results = render(output / "figures", output / "validation")
    if args.check_vectors:
        before = output / "before"
        before.mkdir(exist_ok=True)
        destinations = {
            "oracle_estimated_nuisance_comparison": ROOT / "assets/main",
            "heterogeneous_gap_design_comparison_lines": ROOT
            / "assets/supplement/v3/figures",
            **{
                f"supp_scenario_{scenario}": ROOT / "assets/supplement/v4/figures"
                for scenario in REFERENCES["extended"]
            },
        }
        for stem, folder in destinations.items():
            shutil.copy2(folder / (stem + ".pdf"), before / (stem + ".pdf"))
        subprocess.run(
            [args.pdf_python, str(HERE / "verify_vector_geometry.py"), str(output)],
            check=True,
        )
    inputs_after = {
        p.relative_to(INPUTS).as_posix(): sha(p)
        for p in INPUTS.rglob("*")
        if p.is_file()
    }
    assert inputs_before == inputs_after
    report = {
        "success": True,
        "input_sha256": inputs_before,
        "curve_count": 192,
        "figure_count": 6,
        "vector_check_requested": args.check_vectors,
        "scope": "Figures regenerated from saved summaries.",
    }
    (output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
