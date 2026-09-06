# BAFFLES Changelog

Newest first. The code and grids used in the published paper are at commit `d2ce435` (see the 2020-06-16 entry at the bottom for how to check it out).

## Unreleased: installable Python package (v1.1.0)

Pull request: https://github.com/adamstanfordmoore/BAFFLES/pull/4 (inspired by Luke Bouma's [PR #2](https://github.com/adamstanfordmoore/BAFFLES/pull/2))

### Layout

- All modules moved into a `baffles/` package. `baffles.py` became `baffles/core.py`; `baffles/__init__.py` re-exports `baffles_age`, `age_estimator` and `posterior`, so `import baffles` works exactly as before and the awkward `import baffles.baffles` is avoided.
- Internal imports are absolute (`from baffles import fitting as my_fits`, `from baffles import li_constants as const`, ...) instead of bare top-level module names that polluted `sys.path` and only worked from the repo root.
- `data/` and `grids/` moved inside the package and are shipped as package data. New `baffles/paths.py` exposes `DATA_DIR` and `GRID_DIR`; every hard-coded relative path (`join('data', ...)`, `'grids/...'`) now resolves against them, so the package works from any working directory. The `DEFAULT_MEDIAN_GRID` constants are built from `GRID_DIR`, and `age_estimator.set_default_grids` (used by `refresh.py`) edits the packaged constants file.
- The command-line block from `baffles.py` became `baffles/cli.py` with a `main()`; installed as the `baffles` console script with the same flags (`baffles -bmv 0.65 -rhk -4.906 -plot`).
- Paper and maintenance scripts moved to `scripts/` (`refresh.py`, `paper_plots_li.py`, `paper_plots_ca.py`, `make_baffles_table2.py`, `make_baffles_table3.py`); `usage_example.py` to `examples/`.
- `pyproject.toml` (setuptools) replaces `requirements.txt`; `environment.yml` installs the package in editable mode. Version 1.1.0.
- `tests/test_regression.py`: six pytest checks on reference posteriors (package data present, Sun from calcium, lithium detection, upper limit, `maxAge`, combined posterior).

### Testing performed

- Six reference `baffles_age` cases give byte-identical printed output to the pre-package code.
- `pip install .` into a fresh virtual environment, then running the `baffles` command and `import baffles` from an unrelated directory: all 49 data files and 4 grids are shipped and found.
- `examples/usage_example.py`, `pytest`, and `compileall` with warnings as errors all pass.
- All fifteen manuscript plotting functions run from `scripts/`.
- `scripts/refresh.py` regenerates grids into the package directory and rewrites the constants correctly.

No numerical results changed.

## 2026-09-06: NumPy 2 / SciPy 1.14+ support, lithium residual fix, regenerated grids

Pull request: https://github.com/adamstanfordmoore/BAFFLES/pull/3

### Code changes

Dependencies had been capped at `numpy<2.0` and `scipy<1.13` because newer releases removed functions BAFFLES used. Those calls were replaced and the caps dropped.

| Removed / renamed API | Replacement | Files |
|---|---|---|
| `scipy.interpolate.interp2d` (removed in SciPy 1.14) | `scipy.interpolate.RectBivariateSpline(kx=1, ky=1)` via new helper `fitting.median_grid_interpolator` | `baffles.py`, `fitting.py` |
| `numpy.trapz` (deprecated in NumPy 2.0) | `scipy.integrate.trapezoid` | `baffles.py`, `probability.py` |
| `scipy.integrate.cumtrapz` (removed in SciPy 1.14) | `scipy.integrate.cumulative_trapezoid` | `probability.py`, `fitting.py` |
| `numpy.float` (removed in NumPy 1.24) | `float` | `readData.py`, `paper_plots_ca.py` |
| `matplotlib.cm.get_cmap` (removed in Matplotlib 3.9) | `matplotlib.pyplot.get_cmap` | `plotting.py`, `paper_plots_li.py`, `paper_plots_ca.py` |

The `RectBivariateSpline` with first-order splines is bilinear interpolation on the same grid `interp2d` used, and reproduces `interp2d` output on a grid to 4e-16.

Also: LaTeX backslashes in five label strings escaped (Python 3.12 `SyntaxWarning`), `environment.yml` now targets Python 3.12, `.gitignore` added (`__pycache__/`, `plots/`), README download instructions switched to a blobless clone.

### Bug fix: misaligned lithium residuals

`interp2d.__call__` silently sorted its input coordinates. In `fitting.get_fit_residuals` the per-cluster B-V arrays are in catalog order, not sorted, so each lithium residual was computed as `log(Li EW)_i - i(t, sorted_bv_i)`, pairing every star with a different star's B-V. The replacement evaluates pointwise, keeping each residual aligned with its own star.

- The posterior code path in `baffles.likelihood` was unaffected: its B-V samples from `gaussian_cdf_space` are already ascending, so the sort was a no-op there.
- Calcium is unaffected: its residuals use a one-dimensional interpolation in age only.
- Effect on the lithium residual distribution used to build the likelihood `K(l | i(t,b))`: standard deviation 0.278 dex -> 0.200 dex, median absolute deviation 0.119 -> 0.086 dex. The 2020 lithium likelihood was therefore wider than the data support.

### Regenerated grids

`refresh.py` was rerun on the new code. Files in `grids/` are now dated `090626` and `ca_constants.py` / `li_constants.py` point at them; the `061620` files were removed (retrievable at commit `d2ce435`).

| Grid | Change vs. 2020 |
|---|---|
| `median_rhk_090626.npy` | max abs diff 2e-14 |
| `median_li_090626.npy` | max abs diff 5e-6 |
| `calcium_likelihood_fit.npy` | max abs diff 7e-6 |
| `lithium_likelihood_fit.npy` | rebuilt from aligned residuals; PDF peak 2.49 -> 3.29 |

The cached data pickles in `data/` were rewritten by `refresh.py` but are content-identical (same stars, same fits) and were left as committed.

### Testing performed

Two isolated environments were built with `uv` on Python 3.12:

- Old stack: NumPy 1.26.4, SciPy 1.12.0 (the last versions satisfying the previous caps).
- New stack: NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.1, Astropy 8.0.1 (current releases).

Checks, all run with `DeprecationWarning` and `FutureWarning` promoted to errors:

1. **Code change alone (old stack, shipped 2020 grids).** Six end-to-end `baffles_age` cases (calcium, lithium detection, upper limit, `maxAge`, large B-V error, combined), four direct `likelihood()` evaluations, and the calcium fit residuals were compared between the original code and the new code. Worst relative difference 2e-14. Only the lithium residuals differed, as intended.
2. **Library change (new code, old vs. new stack).** Same cases: worst relative difference 4e-12.
3. **Interpolator equivalence.** `RectBivariateSpline` vs. `interp2d` on the median grid: max abs diff 4e-16. Demonstrated that `interp2d` reorders unsorted inputs and the replacement does not.
4. **`refresh.py` full pipeline** on the new stack, including reading raw data files with Astropy 8: median grids reproduce the 2020 grids to 5e-6.
5. **Final grids on both stacks.** The six `baffles_age` cases give byte-identical printed output on the old and new stacks.
6. **All fifteen manuscript figure functions** in `paper_plots_li.py` and `paper_plots_ca.py` run to completion on the new stack (output in `plots/`, untracked).
7. **`python -m compileall`** with `SyntaxWarning` as error passes.

### Impact on results

Calcium ages are unchanged. Lithium ages shift by a few percent and their credible intervals narrow. Comparison of the regenerated manuscript figures against the published ones (legend values; ranges are 68% intervals):

| Figure | Object | Isochronal / literature age | Paper (2020 grids) | Current |
|---|---|---|---|---|
| 10 | all six calcium clusters | - | - | identical |
| 11 | NGC 2264 | 5 Myr | 4.66 Myr (3.8-5.3) | 4.81 Myr (4.1-5.3) |
| 11 | α Per | 85 Myr | 106 Myr (92-120) | 105 Myr (92-120) |
| 11 | Pleiades | 130 Myr | 124 Myr (110-130) | 125 Myr (120-130) |
| 11 | M35 | 200 Myr | 194 Myr (190-200) | 195 Myr (190-200) |
| 11 | M34 | 240 Myr | 237 Myr (210-250) | 230 Myr (210-250) |
| 11 | Hyades | 700 Myr | 953 Myr (830-1100) | 953 Myr (840-1100) |
| 12 | AB Dor | 149 Myr | 128 Myr (99-160) | 125 Myr (100-150) |
| 12 | Tuc/Hor | 45 Myr | 35.1 Myr (24-47) | 36.3 Myr (27-46) |
| 13 | TW PsA, Li only | 440 ± 40 Myr (lit.) | 295 Myr (210-370) | 297 Myr (240-370) |
| 13 | TW PsA with Fomalhaut PDF | 440 ± 40 Myr (lit.) | 356 Myr (280-410) | 347 Myr (280-390) |
| 13 | HR 2562, combined | 300-900 Myr (lit.) | 658 Myr (520-1100) | 642 Myr (530-1000) |
| 13 | HD 206893, combined | 200-2100 Myr (lit.) | 573 Myr (380-1000) | 567 Myr (390-1000) |

Figure 9 (lithium residual histogram and `K(l|i(t,b))`) is visibly narrower: residuals now lie mostly within ±0.6 dex rather than ±1.2 dex, and the PDF peak rises from about 2.5 to 3.3. Figures 1-8 and 14 depend only on the median grids or on calcium and are unchanged.

Agreement with isochronal ages is unchanged: the same clusters contain the isochronal age within their 68% interval as before, the Hyades remains the outlier, and both moving groups stay within one sigma. Manuscript text values that no longer match the code: AB Dor 127 +35/-28 Myr (now 125 +29/-24), Tuc/Hor 35 +11/-10 Myr (now 36 +9/-8), HD 206893 final age 570 Myr, 380-1000 (now 567, 390-1000).

## 2020-06-16: Published paper version

The code and grids used in [Stanford-Moore et al. 2020](https://arxiv.org/abs/2006.04811) are at commit `d2ce435` (also archived on [Zenodo](https://doi.org/10.5281/zenodo.3840244)). After cloning as described in the README, run:

```bash
git checkout d2ce435
```

That version requires `numpy<2.0` and `scipy<1.13`.
