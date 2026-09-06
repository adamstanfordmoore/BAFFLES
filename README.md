# BAFFLES: Bayesian Ages for Lower-mass Field Stars

This package computes age posteriors for field stars from measurements of R'HK calcium emission and/or B-V color and lithium equivalent width absorption (Li EW). For calcium emission our method is calibrated to stars with B-V between 0.45 and 0.9 (~ F6-K2) and
log(R'HK) between -3.7 and -5. For lithium we have calibrated BAFFLES to stars with B-V between 0.35 and 1.9 (~F2-M5) and Li EW between 3.2 and 1500 mA. See the paper [Stanford-Moore et al. 2020](https://arxiv.org/abs/2006.04811).

### Downloading (Size ~8MB)

Download the zipped file from GitHub or [Zenodo](https://doi.org/10.5281/zenodo.3840244).

Or

`git clone --depth=1 https://github.com/adamstanfordmoore/BAFFLES.git`

Currently there are still some large files in the git history, so clone with depth 1 to avoid all 256MB of history.

### Installation

#### Using Conda (Recommended)

Create a conda environment with all dependencies:

```bash
conda env create -f environment.yml
conda activate baffles
```

#### Using pip

```bash
pip install -r requirements.txt
```

### Requirements

- Python 3.9+ (3.12 recommended)
- numpy >= 1.20 (NumPy 2.x supported)
- scipy >= 1.7
- matplotlib >= 3.3
- astropy >= 4.0

Tested with NumPy 1.26 / SciPy 1.12 and NumPy 2.5 / SciPy 1.18; posteriors agree to machine precision across both stacks.

## Changelog

### 2026-09-06: NumPy 2 / SciPy 1.14+ support and regenerated lithium grids

**Code.** Replaced the removed `scipy.interpolate.interp2d` with a bilinear `RectBivariateSpline`, updated other removed NumPy/SciPy functions (`np.trapz`, `cumtrapz`, `np.float`), and dropped the `numpy<2.0` / `scipy<1.13` caps.

**Bug fix.** `interp2d` silently sorted its input coordinates. In `fitting.get_fit_residuals` the cluster B-V arrays are unsorted, so each lithium residual was paired with the wrong star's B-V, inflating the residual scatter used to build the lithium likelihood (std 0.28 dex; correctly aligned it is 0.20 dex). The posterior code path was unaffected because its B-V samples were already sorted. Calcium is unaffected.

**Impact.** Grids in `grids/` were regenerated with `refresh.py` (now dated `090626`). Median indicator grids are unchanged to within 5e-6. Calcium ages are unchanged. Lithium ages shift modestly and their credible intervals narrow, e.g.:

| Input | Paper version (2020 grids) | Current |
|---|---|---|
| B-V 0.80, Li EW 100 mA | 246 Myr, 68% CI 148-389 | 239 Myr, 68% CI 156-317 |
| B-V 0.45 +/- 0.02, Li EW 21 +/- 5 mA | 774 Myr, 68% CI 499-3060 | 709 Myr, 68% CI 512-2090 |

**Reproducing the published paper results.** The grids and code used in Stanford-Moore et al. 2020 are at commit `d2ce435` (also archived on [Zenodo](https://doi.org/10.5281/zenodo.3840244)). A `--depth=1` clone does not include it, so clone the full history:

```bash
git clone https://github.com/adamstanfordmoore/BAFFLES.git
cd BAFFLES
git checkout d2ce435
```

That version requires `numpy<2.0` and `scipy<1.13`.

## Authors

- Adam Stanford-Moore
- Eric Nielsen
- Rob De Rosa
- Bruce Macintosh
- Ian Czekala

## Command Line Usage Examples

Quick usage from the command line. Let's find an age for the sun using B-V of 0.65 and logR'HK = -4.908 (Mamajek & Hillenbrand 2008).

`python baffles.py -bmv 0.65 -rhk -4.906 -plot`

To save this probability density function in a csv file as 1000 lines of age,probability with optional filename (\_calcium.csv will be appended to name):

`python baffles.py -bmv 0.65 -rhk -4.906 -plot -s -filename suns_age`

Now lets find the age of HR 2562 using B-V=.45 ± .02, log(R'HK) = -4.55 (Gray 2006), and lithium EW of 21 ± 5 (Mesa el al 2018). "-ul" would denote an upper limit. "-s" will save a csv file of the posterior. "-plot" will show a plot of the posterior. -maxAge 10000 will constrain the prior on age to be uniform out to 10 Gyr. "-li_err" allows input of uncertainty
on Li EW, and "-bv_err" uncertainty on B-V. The following command will determine the age using calcium and lithium separately and then find the combined posterior product.

`python baffles.py -bmv 0.45 -bmv_err .02 -rhk -4.55 -li 21 -li_err 5 -plot`

Type `python baffles.py -help` into the command line to learn more.

To directly import baffles and use the module in a python script see "usage_example.py"

## Package Files

**baffles.py** : main file that computes age posteriors. It makes and stores grids of mean indicator as functions of age and B-V

**usage_example.py** : example script using the baffles package

**refresh.py** : updates any stored fits or grids (like mean R'HK as a function of age). Run `python refresh.py` whenever constants or data are changed.

**ca_constants.py** : constants related to calcium

**li_constants.py** : constants related to lithium

**fitting.py** : various fitting functions used to compute grids in baffles.py, inlcuding for mean R'HK as a function of age and mean LiEW as a function of age and B-V

**plotting.py** : plotting functions to display posteriors, data, and fits

**probability.py** : assortment of statistical functions like finding gaussian PDF/CDF

**readData.py** : reads in lithium and calcium data from data directory and returns data and indicator vs B-V fits

**utils.py** : extra helper functions

**paper_plots_li.py** : more advanced plotting examples for lithium used to create the plots in stanford_moore et al 2019.

**paper_plots_ca.py** : more advanced plotting for calcium

**data/** : directory contains files with calcium/lithium data and .p pickle files with saved arrays of indicator vs B-V

**grids/** : directory contains saved grids of mean/sigma indicator values as functions of age/B-V. Refresh these grids with refresh.py
