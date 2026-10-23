# Anytime-valid weighted survival tests for exact failures and completed assessment windows

Code, data and numerical results for the paper by Hien Duy Nguyen and Liming Xiang.

## Manuscript version

The exact 62-page revised manuscript and its LaTeX source are archived at
[`manuscript-round7-2026-10-23`](https://github.com/hiendn/Anytime_Valid_Survival_Tests/tree/manuscript-round7-2026-10-23).
The scientific analysis code, data and saved numerical results remain those at
[`025903c484a004cd278c86738f0721d05104c790`](https://github.com/hiendn/Anytime_Valid_Survival_Tests/tree/025903c484a004cd278c86738f0721d05104c790).
The manuscript archive is a later version: it adds the independent-stream
sufficient result, clarifies complete-specification and pilot-conditional
comparisons, updates literature and changes the monitoring figure's labeling.
No simulation results or patient records were changed or regenerated.

| Archived manuscript file | SHA-256 |
| --- | --- |
| `main.tex` | `6110c5a1c9a3c2b8bd2590148e98477367c901b9d721a1d5e59d20660cefc840` |
| `manuscript.pdf` | `d521bc2fca7fe55a813ae8d2272b9fa2794c8d9b9e7124722d6b64f6490bce53` |

`MANIFEST.json` records this manuscript version and the source/PDF target hashes.
`Verify_Results.py` checks those hashes before checking tables and saved results;
its output identifies the exact source and PDF checked. `VERIFICATION.json`
records the successful run for this archive. A replay of saved comparisons is
not a new simulation run or validation of the clinical observation assumptions.

## Requirements

The reference studies used R 4.4.0. The programs use `data.table`, `jsonlite`, `Rcpp`, `matrixStats`, `digest`, `survival`, `interval`, `msm` and `ggplot2`. The comparator tests also use `testthat` and a C++17 compiler. Install the R packages with:

```r
install.packages(c("data.table", "jsonlite", "Rcpp", "matrixStats", "digest",
                  "survival", "interval", "msm", "ggplot2", "testthat"))
```

Python 3 is sufficient to verify the saved results. Building the manuscript requires LaTeX.

## Reproduce and verify

Run from the repository root:

```sh
python3 Verify_Results.py
Rscript --vanilla Reproduce_Figures.R
```

The first command independently recalculates all 1,140 saved clinical evidence values from the comparison tables, checks the corresponding manuscript summaries and checks the fixed-design null results. The second recreates the four manuscript figures from the saved numerical results. Neither command reruns the simulations or patient preprocessing.

The optional comparator check is:

```sh
Rscript --vanilla tests/session/test-sim2.R
```

It checks the fixed and prequential calculations on deterministic tables, including prediction before simultaneous releases and null normalization. It compiles the numerical kernel and stops if a test fails.

## Source files

| File or directory | Purpose |
| --- | --- |
| `Verify_Results.py` | Independent checks of saved clinical paths and manuscript summaries |
| `Reproduce_Figures.R` | Recreate the four manuscript figures |
| `Clinical_Figures.R` | Plot clinical evidence from saved paths |
| `Export_Data.R` | Export the public clinical example data from `msm` |
| `simulations/code/` | Event generation, comparison tables, simulations and clinical construction |
| `scripts/run_remaining_simulations.R` | Fixed specifications, mixtures, prequential procedures and interval-censored benchmark |
| `scripts/sim2_comparators.R` and `scripts/sim2_kernels.cpp` | Conditional likelihoods and prequential fitting |
| `tests/session/` | Comparator implementation checks |

## Reference outputs

`simulations/results/` contains the primary, weighting, observation-scheme and sensitivity results. The fixed-specification results, paired Monte Carlo uncertainties, full-precision settings and seed records are in `results/session/todo_production_2026-10-15/`. Clinical comparison tables and evidence paths are in `reference/liming_revision/data/`; the diagnostics and visit-time stress results are in `results/session/referee_2026-10-22/clinical/`.

The simulation programs use `L'Ecuyer-CMRG` random streams. The broad simulation battery uses seed `20260901`; the fixed-specification study records separate null-pilot and evaluation seeds. Exact and interval analyses of a simulated dataset share its latent event and censoring times. `MANIFEST.json` records the committed file checksums.

Full simulation reruns are separate from the checks above and take substantially longer. Use fresh output directories:

```sh
Rscript --vanilla simulations/code/run_practical_battery.R --output-dir simulations/results_new
Rscript --vanilla scripts/run_remaining_simulations.R --mode production --output results/session/new_fixed_run
```

The additional covariate, mixed-observation and clinical drivers are in `simulations/code/`. Their source gives the design settings and command-line options. The clinical driver can write a fresh analysis without changing the saved reference data:

```sh
Rscript --vanilla simulations/code/run_real_data_applications.R --output-dir results/clinical_new --skip-figures
```

## Data

The clinical examples use the public cardiac allograft vasculopathy and psoriatic arthritis data distributed with the R package [`msm`](https://chjackson.github.io/msm/). CSV copies from `msm` 1.8.2 are included in `data/`, with their source and licence. `Rscript --vanilla Export_Data.R` recreates these exports. The clinical programs use the package data directly.

The reconstructed comparisons and clinical threshold crossings rely on the working observation assumptions described in the paper. Recalculating evidence from the exported comparison tables does not validate those assumptions or independently validate patient preprocessing.

## Contact

Hien Duy Nguyen: h.nguyen5@latrobe.edu.au
