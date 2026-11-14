"""Reproduce the paper's eight planning figures and their numeric summaries."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--rscript", help="Rscript executable; otherwise discovered on PATH"
    )
    parser.add_argument(
        "--skip-r", action="store_true", help="omit the sampling-distribution plot saved-draw replot"
    )
    parser.add_argument(
        "--font-family",
        default="sans",
        help="R graphics font for the sampling-distribution plot (default: sans; paper font: Arial)",
    )
    args = parser.parse_args()
    rscript = None
    if not args.skip_r:
        rscript = args.rscript or shutil.which("Rscript")
        if not rscript:
            parser.error(
                "Rscript is needed for the sampling-distribution plot; supply --rscript or use --skip-r."
            )
        if args.rscript and not (
            Path(args.rscript).is_file() or shutil.which(args.rscript)
        ):
            parser.error(f"Rscript executable not found: {args.rscript}")
    if args.output_dir:
        os.environ["HCT_PLANNING_OUTPUT"] = str(args.output_dir.resolve())
    from planning_paths import output_root

    output = output_root()
    source = Path(__file__).resolve().parent
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "MPLBACKEND": "Agg"}
    for script in (
        "plot_zero_gap_power_publication.py",
        "plot_study_design_sensitivity.py",
        "common_discrepancy.py",
        "gap_analytical_figures.py",
    ):
        subprocess.run([sys.executable, str(source / script)], env=env, check=True)
    if not args.skip_r:
        subprocess.run(
            [
                rscript,
                str(source / "illustrate_hct_targets_and_sampling.R"),
                "--output-dir",
                str(output),
                "--font-family",
                args.font_family,
            ],
            env=env,
            check=True,
        )
    print(f"Planning outputs: {output}")


if __name__ == "__main__":
    main()
