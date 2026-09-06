"""Regression tests: posterior statistics for reference stars.

Expected values were produced with the grids dated 090626 (NumPy 2.5 / SciPy 1.18)
and verified identical on NumPy 1.26 / SciPy 1.12.  Stats are ages (Myr) at the
CDF quantiles [0.023, 0.16, 0.5, 0.84, 0.977]; for upper limits [0.0027, 0.046, 0.32, 1].
"""
import os
import numpy as np
import pytest
from scipy.integrate import trapezoid

import baffles
from baffles import ca_constants, li_constants


def test_package_data_present():
    assert os.path.isfile(ca_constants.DEFAULT_MEDIAN_GRID)
    assert os.path.isfile(li_constants.DEFAULT_MEDIAN_GRID)
    assert np.load(li_constants.DEFAULT_MEDIAN_GRID).shape == (len(li_constants.BV_S), len(li_constants.AGE))


def test_sun_calcium():
    p = baffles.age_estimator('calcium').get_posterior(0.65, -4.906)
    assert p.stats[2] == pytest.approx(7080, rel=0.01)
    assert p.stats[1] == pytest.approx(4230, rel=0.01)
    assert p.stats[3] == pytest.approx(10300, rel=0.01)
    assert p.array.shape == ca_constants.AGE.shape
    assert trapezoid(p.array, ca_constants.AGE) == pytest.approx(1.0, abs=1e-6)


def test_lithium_detection():
    p = baffles.age_estimator('lithium').get_posterior(0.8, 100)
    assert p.stats[2] == pytest.approx(239, rel=0.01)
    assert p.stats[1] == pytest.approx(156, rel=0.01)
    assert p.stats[3] == pytest.approx(317, rel=0.01)


def test_lithium_with_uncertainties_and_max_age():
    p = baffles.age_estimator('lithium').get_posterior(0.55, 50, bv_uncertainty=None,
                                                      measure_err=None, upperLim=False, maxAge=1000)
    assert p.stats[2] == pytest.approx(665, rel=0.01)
    assert p.stats[4] <= 1000 + 1e-6


def test_lithium_upper_limit():
    p = baffles.age_estimator('lithium').get_posterior(0.7, 10, upperLim=True)
    assert p.upperLim
    assert p.stats[2] == pytest.approx(6140, rel=0.01)  # 1-sigma lower limit
    assert p.stats[1] == pytest.approx(2860, rel=0.01)  # 2-sigma lower limit


def test_baffles_age_combined(capsys):
    p = baffles.baffles_age(bv=0.45, rhk=-4.55, li=21, bv_err=.02, li_err=5, showPlots=False)
    out = capsys.readouterr().out
    assert "Final Median Age: 686 Myr" in out
    assert p.stats[2] == pytest.approx(686, rel=0.01)
