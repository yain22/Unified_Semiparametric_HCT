# Simulation figures

From the repository root:

```sh
python analysis/simulation_figures/plot_simulations.py
```

This writes six PDFs to `results/simulation-figures/figures/` from the saved summaries in `inputs/`. It checks all 192 curve-coordinate hashes against the reference results. It does not rerun any simulations.

`--output-dir PATH` selects another destination. `--check-vectors` additionally compares the plot geometry with the reference PDFs in `assets/` and requires PyMuPDF. Use `--pdf-python PATH` if PyMuPDF is installed in another Python environment.

For newly generated 5,000-repetition results, use `analysis/simulations/plot_publication.py`. The publication renderers require complete summaries, so they should not be applied to four-repetition smoke runs.

