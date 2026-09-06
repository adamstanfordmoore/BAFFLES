# BAFFLES: Bayesian Ages for Lower-mass Field Stars

This package computes age posteriors for field stars from measurements of R'HK calcium emission and/or B-V color and lithium equivalent width absorption (Li EW). For calcium emission our method is calibrated to stars with B-V between 0.45 and 0.9 (~ F6-K2) and
log(R'HK) between -3.7 and -5. For lithium we have calibrated BAFFLES to stars with B-V between 0.35 and 1.9 (~F2-M5) and Li EW between 3.2 and 1500 mA. See the paper [Stanford-Moore et al. 2020](https://arxiv.org/abs/2006.04811).

## Downloading (Size ~12MB)

```bash
git clone --filter=blob:none https://github.com/adamstanfordmoore/BAFFLES.git
```

The `--filter=blob:none` option is recommended because the git history contains large grid files (about 770MB) that were deleted long ago. A blobless clone downloads all current files and the full commit history, but skips old file contents unless you check out an old commit, so the download is only a few MB. A plain `git clone` also works but downloads the whole history.

Alternatively, download the zipped file from GitHub or [Zenodo](https://doi.org/10.5281/zenodo.3840244).

## Installation

BAFFLES is a Python package. Install it into your environment from the cloned directory:

```bash
pip install .
```

Use `pip install -e .` for an editable install if you plan to modify the code or regenerate grids. The calibration data and grids are shipped inside the package, so it works from any directory once installed.

### Using Conda

```bash
conda env create -f environment.yml
conda activate baffles
```

This creates the environment and installs BAFFLES into it in editable mode.

## Requirements

- Python 3.9+ (3.12 recommended)
- numpy >= 1.20 (NumPy 2.x supported)
- scipy >= 1.7
- matplotlib >= 3.3
- astropy >= 4.0

Tested with NumPy 1.26 / SciPy 1.12 and NumPy 2.5 / SciPy 1.18; posteriors agree to machine precision across both stacks. Run the regression tests with `pip install pytest && pytest`.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for details of changes since publication, the testing performed, and their impact on derived ages. The code and grids used in the published paper are at commit `d2ce435`; after cloning as described above, run `git checkout d2ce435` to view them.

## Authors

- Adam Stanford-Moore
- Eric Nielsen
- Rob De Rosa
- Bruce Macintosh
- Ian Czekala

## Command Line Usage Examples

Installing the package provides a `baffles` command. Let's find an age for the sun using B-V of 0.65 and logR'HK = -4.906 (Mamajek & Hillenbrand 2008).

`baffles -bmv 0.65 -rhk -4.906 -plot`

To save this probability density function in a csv file as 1000 lines of age,probability with optional filename (\_calcium.csv will be appended to name):

`baffles -bmv 0.65 -rhk -4.906 -plot -s -filename suns_age`

Now let's find the age of HR 2562 using B-V=.45 ± .02, log(R'HK) = -4.55 (Gray 2006), and lithium EW of 21 ± 5 (Mesa et al. 2018). "-ul" would denote an upper limit. "-s" will save a csv file of the posterior. "-plot" will show a plot of the posterior. -maxAge 10000 will constrain the prior on age to be uniform out to 10 Gyr. "-li_err" allows input of uncertainty
on Li EW, and "-bv_err" uncertainty on B-V. The following command will determine the age using calcium and lithium separately and then find the combined posterior product.

`baffles -bmv 0.45 -bmv_err .02 -rhk -4.55 -li 21 -li_err 5 -plot`

Type `baffles -help` into the command line to learn more.

## Python Usage

```python
import baffles

# one line, like the command line
posterior = baffles.baffles_age(bv=0.65, rhk=-4.906)

# or build a reusable estimator
est = baffles.age_estimator('lithium')
posterior = est.get_posterior(bv=0.8, metallicity=100)   # Li EW in mA
posterior.array   # PDF on the age grid baffles.li_constants.AGE
posterior.stats   # ages at CDF [.02, .16, .5, .84, .97]
```

See `examples/usage_example.py` for posterior products and plotting. Lower-level functionality is in the submodules `baffles.fitting`, `baffles.plotting`, `baffles.probability`, `baffles.readData`, `baffles.utils`, `baffles.ca_constants` and `baffles.li_constants`.

## Repository Layout

**baffles/** : the installable package

- **core.py** : computes age posteriors (`baffles_age`, `age_estimator`, `posterior`); makes and stores grids of median indicator as functions of age and B-V
- **cli.py** : the `baffles` command-line entry point
- **ca_constants.py**, **li_constants.py** : constants related to calcium and lithium, including which grid files are the defaults
- **fitting.py** : fitting functions used to compute grids, including mean R'HK as a function of age and mean Li EW as a function of age and B-V
- **plotting.py** : plotting functions to display posteriors, data, and fits
- **probability.py** : statistical helpers such as gaussian PDF/CDF
- **readData.py** : reads the lithium and calcium calibration data and returns data and indicator vs B-V fits
- **utils.py** : extra helper functions
- **paths.py** : locations of the packaged `data/` and `grids/` directories
- **data/** : calcium/lithium calibration data and .p pickle files with cached arrays of indicator vs B-V
- **grids/** : saved grids of median indicator values as functions of age/B-V, and the fitted likelihood functions

**scripts/** : maintenance and paper scripts, run with `python scripts/<name>.py` after installing the package

- **refresh.py** : regenerates all stored fits and grids inside the package. Run it whenever constants or data are changed.
- **paper_plots_li.py**, **paper_plots_ca.py** : produce the figures in Stanford-Moore et al. 2020 (written to `plots/`)
- **make_baffles_table2.py**, **make_baffles_table3.py** : produce the paper's LaTeX tables

**examples/usage_example.py** : example script using the package

**tests/** : pytest regression tests on reference posteriors
