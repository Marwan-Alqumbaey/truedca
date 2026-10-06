# truedca (Python)

Decision curve analysis for outcomes recorded with error. This is the Python
version of the R package in [`../r`](../r); the tests check that both give the
same numbers.

```sh
pip install "git+https://github.com/Marwan-Alqumbaey/truedca#subdirectory=python"
```

Requires Python 3.9 or later with NumPy, SciPy and pandas (`requirements.txt`);
`plot_dca` also needs Matplotlib.

Functions: `nb_curve`, `dca_window`, `warp_threshold`, `nb_adjust`, `nb_bounds`,
`nb_bias_analysis`, `pilot_v`, `twophase_design`, `design_gain`, `twophase_dca`,
`make_strata`, `plot_dca`. They take NumPy arrays and return pandas data frames.
Functions that use random numbers take `rng` (a seed or `numpy.random.Generator`).
In `twophase_dca`, `truth` is `NaN` for patients whose chart was not reviewed.

## Tests

```sh
pip install -r requirements-test.txt
pip install .
pytest
```

The reference values in `tests/fixtures` come from the R package
(`Rscript tests/fixtures/make_fixtures.R`, run from this folder).
