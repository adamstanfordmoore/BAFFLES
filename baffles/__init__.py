"""
BAFFLES: Bayesian Ages for Field LowEr-mass Stars

Stanford-Moore, Nielsen, De Rosa, Macintosh & Czekala (2020), ApJ.
https://arxiv.org/abs/2006.04811

Quick start::

    import baffles
    post = baffles.baffles_age(bv=0.65, rhk=-4.906)   # the Sun from R'HK
    est = baffles.age_estimator('lithium')             # reusable estimator
    post = est.get_posterior(bv=0.8, metallicity=100)  # Li EW in mA

Submodules: baffles.fitting, baffles.plotting, baffles.probability,
baffles.readData, baffles.utils, baffles.ca_constants, baffles.li_constants.
"""
from baffles.core import baffles_age, age_estimator, posterior
from baffles.paths import DATA_DIR, GRID_DIR

__version__ = "1.1.0"
__all__ = ["baffles_age", "age_estimator", "posterior", "DATA_DIR", "GRID_DIR", "__version__"]
