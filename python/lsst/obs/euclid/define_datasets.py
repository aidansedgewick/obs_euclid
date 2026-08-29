from lsst.daf.butler import Butler, DatasetType

from lsst.obs.euclid import euclid_constants

STANDARD_DIMENSIONS = ("skymap", "tract", "patch", "band")


def register_dataset_types(butler: Butler) -> None:
    universe = butler.dimensions

    # Define and register BGMOD (BackGround MODel)
    bgmod_dataset_type = DatasetType(
        euclid_constants.EUCLID_BGMOD_DATASET_TYPE_NAME,
        dimensions=STANDARD_DIMENSIONS,
        storageClass="ExposureF",
        universe=universe,
    )
    butler.registry.registerDatasetType(bgmod_dataset_type)

    catalog_psf_dataset_type = DatasetType(
        euclid_constants.EUCLID_CATALOGPSF_DATASET_TYPE_NAME,
        dimensions=STANDARD_DIMENSIONS,
        storageClass="AstropyTable",
        universe=universe,
    )
    butler.registry.registerDatasetType(catalog_psf_dataset_type)

    # Define and register COADDS
    coadd_dataset_type = DatasetType(
        euclid_constants.EUCLID_COADD_DATASET_TYPE_NAME,
        dimensions=STANDARD_DIMENSIONS,
        storageClass="ExposureF",
        universe=universe,
    )
    butler.registry.registerDatasetType(coadd_dataset_type)

    # Define and register FLAGMAP
    flagmap_dataset_type = DatasetType(
        euclid_constants.EUCLID_FLAGMAP_DATASET_TYPE_NAME,
        dimensions=STANDARD_DIMENSIONS,
        storageClass="ExposureI",  # ExposureI as flags are int32!
        universe=universe,
    )
    butler.registry.registerDatasetType(flagmap_dataset_type)

    # Define and register GRIDPSF
    gridpsf_dataset_type = DatasetType(
        euclid_constants.EUCLID_GRIDPSF_DATASET_TYPE_NAME,
        dimensions=STANDARD_DIMENSIONS,
        storageClass="ExposureF",
        universe=universe,
    )
    butler.registry.registerDatasetType(gridpsf_dataset_type)

    # Define and register RMSMAP
    rmsmap_dataset_type = DatasetType(
        euclid_constants.EUCLID_RMSMAP_DATASET_TYPE_NAME,
        dimensions=STANDARD_DIMENSIONS,
        storageClass="ExposureF",
        universe=universe,
    )
    butler.registry.registerDatasetType(rmsmap_dataset_type)

    # Define and register SEGMAP
    segmap_dataset_type = DatasetType(
        euclid_constants.EUCLID_SEGMAP_DATASET_TYPE_NAME,
        dimensions=["skymap", "tract", "patch"],  # NO band
        storageClass="NumpyArray",  # they are int64 - can't symlink, duplicate data.
        universe=universe,
    )
    butler.registry.registerDatasetType(segmap_dataset_type)
