import re
from logging import getLogger
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from lsst.daf.butler import Butler
from lsst.skymap import DiscreteSkyMap, DiscreteSkyMapConfig

from lsst.obs.euclid import euclid_constants

logger = getLogger("euclid_skymap")


def build_euclid_q1_skymap_config_from_tiles(
    filelist: list[Path],
    hdu: int,
):

    tract_data = {}

    filelist = sorted(filelist)

    # Capturing groups tile_id and rng_id: ...TILE<tile_id>-<rng_id>_...
    pattern = re.compile(".*_TILE([0-9]*)-([a-zA-Z0-9]*)_.*\.fits")

    for ii, filepath in enumerate(filelist):

        matches = pattern.search(filepath.name)

        tile_id = matches.group(1)
        rng_id = matches.group(2)  # May be useful to extract in future.

        if tile_id in tract_data:
            continue  # Don't care if we already know about this tile.

        with fits.open(filepath) as f:
            header = f[hdu].header
            shape = f[hdu].data.shape
            if shape != (19200, 19200):
                logger.info(f"tile {tile_id} has shape {shape} (tract={len(tract_data)})")

            wcs = WCS(header)
            center_pix = (header["NAXIS1"] / 2.0, header["NAXIS2"] / 2.0)
            center_coord: SkyCoord = wcs.pixel_to_world(*center_pix)

        tract_data[tile_id] = {"ra": center_coord.ra.deg, "dec": center_coord.dec.deg}

    logger.info(f"Build skymap with {len(tract_data)} tracts")

    raList = []  # match LSST naming scheme here...
    decList = []
    for tile_id in sorted(tract_data.keys()):
        raList.append(tract_data[tile_id]["ra"])
        decList.append(tract_data[tile_id]["dec"])

    config = DiscreteSkyMapConfig()
    config.raList = raList
    config.decList = decList
    config.radiusList = [0.25] * len(tract_data)  # R per tile
    config.projection = "TAN"
    config.pixelScale = 0.1
    config.tractOverlap = 0.2 / 60.0

    # Choose tract size as the largest euclid tile.
    # Butler is flexible enough to handle this.
    # tile 102041659 is shape (21600, 19200)
    # tile 102160612 is shape (19200, 21600)
    config.tractBuilder["legacy"].patchInnerDimensions = (21600, 21600)
    config.tractBuilder["legacy"].patchBorder = 100
    return config


def read_skymap_from_butler(butler: Butler):
    return butler.get(
        "skyMap",
        dataId={"skymap": euclid_constants.EUCLID_SKYMAP_NAME},
        collections="skymaps",
    )
