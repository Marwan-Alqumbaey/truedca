# Reproducing the paper

Code and results for: Alqumbaey M, Nasher YM. *Decision curve analysis with misclassified
outcomes: bias, sensitivity analysis and optimal two-phase validation.*

Run every command **from this `paper/` folder**. First install the R packages
(R 4.1 or later) and, for the numerical checks, the Python packages:

```sh
Rscript install.R
pip install -r requirements.txt
```

| Step | Command | Time |
| --- | --- | --- |
| 1. Simulation, part A | `Rscript simulation/run_simulation.R A` | about 1 minute |
| 2. Simulation, part B (the two scenarios can run in parallel) | `Rscript simulation/run_simulation.R M3` and `Rscript simulation/run_simulation.R M4` | about 2 hours each |
| 3. Simulation performance measures | `Rscript simulation/run_simulation.R summary` | seconds |
| 4. Supplementary simulations | `Rscript simulation/run_sensitivity.R scenarios` and `Rscript simulation/run_sensitivity.R design` | about 10 minutes each |
| 5. Example 1, eICU | `Rscript analysis/eicu_aki.R ../data/eicu-demo` | about 3 minutes |
| 6. Examples 2 and 3, NHANES | `Rscript ../data/nhanes/download_nhanes.R` (run from the repository root), then `Rscript analysis/nhanes.R ../data/nhanes` | about 5 minutes |
| 7. Figures 1-4 | `Rscript analysis/figures.R` | seconds |
| 8. Tables 1-3 and S2-S8 | `Rscript tables/make_tables.R` | seconds |
| 9. Numerical checks of the theorems | `python verify/check_theorems_1_2_corollary_1.py` (and the other scripts in `verify/`) | about 5 minutes in total |

Times are for one core of a standard laptop. Every script sets its own random seed, so
each step gives the same numbers on every run. The stored results (`simulation/results/`,
`analysis/results/`, `tables/`, `figures/`) are the ones reported in the paper, so steps 7
and 8 can be run without repeating the slow steps. To check a fresh run against the stored
results, copy a results folder before running and compare, for example:

```sh
cp -r analysis/results expected
Rscript analysis/eicu_aki.R ../data/eicu-demo
Rscript check_results.R expected analysis/results
```

| Folder | Contents |
| --- | --- |
| `simulation/` | Simulation study (ADEMP) and its results |
| `analysis/` | The three examples, the figure script and the results of the examples |
| `tables/` | Tables of the paper and its supplement, built from the result files |
| `figures/` | Figures 1-4 (vector PDF and 800 dpi TIFF) |
| `verify/` | Numerical checks of every theorem and of the split score interval |

| Result in the paper | Script | Result files |
| --- | --- | --- |
| Figure 1 | `simulation/run_simulation.R A` | `simulation/results/partA_distortion.csv` |
| Table 1, Figure 2, Tables S2-S5 | `simulation/run_simulation.R M3`, `M4`, `summary` | `simulation/results/partB_*.csv`, `design_gain_predicted.csv` |
| Table S6 | `simulation/run_sensitivity.R` | `simulation/results/sens_*.csv` |
| Section 7, Figure 3, Tables 2, 3 and S7 | `analysis/eicu_aki.R` | `analysis/results/eicu_*.csv` |
| Section 8, Figure 4, Tables 2, 3 and S8 | `analysis/nhanes.R` | `analysis/results/nhanes_*.csv` |
