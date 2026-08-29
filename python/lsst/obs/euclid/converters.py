from pathlib import Path

import numpy as np

from astropy.io import fits
from astropy.table import Table


def catalog_psf_converter(filepath: Path):
    return Table.read(filepath)

def segmap_converter(filepath: Path):
    with fits.open(filepath) as f:
        segmap = np.asarray(f[0].data, dtype=np.int64)
    return segmap
