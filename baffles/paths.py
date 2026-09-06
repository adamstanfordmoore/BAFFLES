"""Locations of the data and grid files shipped inside the baffles package."""
from os.path import dirname, abspath, join

PACKAGE_DIR = dirname(abspath(__file__))
DATA_DIR = join(PACKAGE_DIR, 'data')
GRID_DIR = join(PACKAGE_DIR, 'grids')
