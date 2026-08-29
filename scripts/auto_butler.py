import json
import logging
import os
import re
import sys
from argparse import ArgumentParser
from datetime import datetime, UTC
from logging.handlers import TimedRotatingFileHandler
from pathlib import Path
from typing import Callable

import tqdm

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from lsst.daf.butler import (
    Butler,
    ButlerConfig,
    CollectionType,
    Config,
    DataCoordinate,
    DatasetRef,
    FileDataset,
)
from lsst.geom import Box2D, SpherePoint, degrees
from lsst.skymap import (
    DiscreteSkyMap,
    DiscreteSkyMapConfig,
    Index2D,
    PatchInfo,
    TractInfo,
)

from lsst.obs.euclid import euclid_constants
from lsst.obs.euclid import converters
from lsst.obs.euclid.define_datasets import register_dataset_types
from lsst.obs.euclid.instrument import EuclidInstrument
from lsst.obs.euclid.skymap import (
    build_euclid_q1_skymap_config_from_tiles,
    read_skymap_from_butler,
)

logger = logging.getLogger("auto_butler")

###===== Some definitions about the structure of the euclid files =====###

ALLOWED_DATASETS = ("catalog-psf", "coadd", "grid-psf", "bgmod", "rms", "flag", "segmap")
ALLOWED_CAMERA_SUBDIRS = ("DECAM", "NISP", "VIS")

MANIFEST_TABLES_PATH = Path(__file__).parent.parent / "ingest_manifests" 
# should get <top-level obs_euclid>/ingest_manifests

DATASET_PARENT_PATH_LOOKUP = {
    "catalog-psf": "MER",
    "coadd": "MER",
    "grid-psf": "MER",
    "bgmod": "MER",
    "rms": "MER",
    "flag": "MER",
    "segmap": "MER_SEG",  # The only one that is not MER - also has no subdirs.
}

# Asterisk in this prefix will 'glob' all of the different bandpasses AND tile names
# eg. EUC_MER_CATALOG-PSF-{NIR-H, NIR-J, NIR-Y}_TILE102020531...
# or  EUC_MER_CATALOG-PSF-VIS_TILE102020531...
# or  EUC_MER_CATALOG-PSF-DES-G_TILE102020531...
GLOB_PATTERN_LOOKUP = {
    "catalog-psf": "EUC_MER_CATALOG-PSF-*.fits",  # one *: for both bpass AND tile_id
    "coadd": "EUC_MER_BGSUB-MOSAIC-*.fits",
    "grid-psf": "EUC_MER_GRID-PSF-*.fits",
    "bgmod": "EUC_MER_BGMOD-*.fits",
    "rms": "EUC_MER_MOSAIC-*-RMS_*.fits",  # two asterisk: 1st bpass, 2nd tile_id
    "flag": "EUC_MER_MOSAIC-*-FLAG_*.fits",  # two asterisks.
    "segmap": "EUC_MER_FINAL-SEGMAP*.fits",
}


# Capturing groups tile_id and rng_id: ...TILE<tile_id>-<rng_id>_...
# grid-psf: '\-': capture literal '-' in eg. NIR-Y
BANDPASS_REGEX_LOOKUP = {
    "grid-psf": r"EUC_MER_GRID-PSF-([A-Z\-]*)_TILE.*",
    "catalog-psf": r"EUC_MER_CATALOG-PSF-([A-Z\-]*)_TILE.*"
}

# The relevant header is different for some datasets
HDU_LOOKUP = {
    **{dataset: 0 for dataset in ["coadd", "bgmod", "rms", "flag", "segmap"]},
    **{dataset: 1 for dataset in ["grid-psf", "catalog-psf"]},
}

HAS_CAMERA_SUBDIRS = ("catalog-psf", "coadd", "grid-psf", "bgmod", "rms", "flag")
# "segmap" does NOT have subdirs.

PHYSICAL_FILTER_TO_BAND = {"VIS": "VIS", "NIR_Y": "Y", "NIR_J": "J", "NIR_H": "H"}

PUT_CONVERTER_LOOKUP = {
    "segmap": converters.segmap_converter,
    "catalog-psf": converters.catalog_psf_converter,
}


dataset_type_lookup = {
    "bgmod": euclid_constants.EUCLID_BGMOD_DATASET_TYPE_NAME,
    "catalog-psf": euclid_constants.EUCLID_CATALOGPSF_DATASET_TYPE_NAME,
    "coadd": euclid_constants.EUCLID_COADD_DATASET_TYPE_NAME,
    "flag": euclid_constants.EUCLID_FLAGMAP_DATASET_TYPE_NAME,
    "grid-psf": euclid_constants.EUCLID_GRIDPSF_DATASET_TYPE_NAME,
    "rms": euclid_constants.EUCLID_RMSMAP_DATASET_TYPE_NAME,
    "segmap": euclid_constants.EUCLID_SEGMAP_DATASET_TYPE_NAME,
}

collection_name_lookup = {
    "bgmod": euclid_constants.EUCLID_BGMOD_COLLECTION_NAME,
    "catalog-psf": euclid_constants.EUCLID_CATALOGPSF_DATASET_TYPE_NAME,
    "coadd": euclid_constants.EUCLID_COADD_COLLECTION_NAME,
    "flag": euclid_constants.EUCLID_FLAGMAP_COLLECTION_NAME,
    "grid-psf": euclid_constants.EUCLID_GRIDPSF_COLLECTION_NAME,
    "rms": euclid_constants.EUCLID_RMSMAP_COLLECTION_NAME,
    "segmap": euclid_constants.EUCLID_SEGMAP_COLLECTION_NAME,
}


###===== define functions here =====###


def collect_files(dataset: str, euclid_data_base_path: Path, cameras=("VIS", "NISP")):

    if dataset not in ALLOWED_DATASETS:
        raise KeyError(f"Unknown dataset '{dataset}': choose from {ALLOWED_DATASETS}")

    parent_path = DATASET_PARENT_PATH_LOOKUP[dataset]
    data_path = euclid_data_base_path / parent_path

    tile_directories = list(sorted(data_path.glob("*")))
    glob_pattern = GLOB_PATTERN_LOOKUP[dataset]

    filelist = []
    for ii, tile_dir in enumerate(tile_directories):
        if dataset in HAS_CAMERA_SUBDIRS:
            for camera in cameras:
                camera_subdir = tile_dir / camera
                files = list(camera_subdir.glob(glob_pattern))
                filelist.extend(files)
        else:
            files = list(tile_dir.glob(glob_pattern))
            filelist.extend(files)

    logger.info(f"collected {len(filelist)} '{dataset}' files")
    return sorted(filelist)


def build_euclid_skymap(
    repo_path: Path,
    euclid_data_base_path: Path,
    tracts_from_dataset: str = "coadd",
    cameras: tuple[str] = ("VIS",),  # Choose VIS only - don't need NISP/DECAM,
    skymap_config_path: Path = None,
):
    dataset_filelist = collect_files(tracts_from_dataset, euclid_data_base_path, cameras=cameras)
    hdu = HDU_LOOKUP[tracts_from_dataset]
    glob_pattern = GLOB_PATTERN_LOOKUP[tracts_from_dataset]

    ##=== Build config from file
    skymap_config: DiscreteSkyMapConfig = build_euclid_q1_skymap_config_from_tiles(
        dataset_filelist, hdu
    )

    ##=== Write the config for inspection later
    skymap_config_path = repo_path / euclid_constants.EUCLID_SKYMAP_FILENAME
    tstr = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S")
    skymap_comments = (
        f"# DiscreteSkyMap of {len(skymap_config.raList)} tracts.\n"
        f"# Tracts are defined per-tile from all tiles in Euclid Q1\n"
        f"#     from all (unique) tiles from files {glob_pattern}.\n"
        f"# Generated at {tstr}UT\n\n"
    )
    skymap_script: str = skymap_config.saveToString(skymap_config_path)

    with open(skymap_config_path, "w+") as f:
        lines = skymap_comments + skymap_script
        f.write(lines)

    ##=== The actual building of the skymap from config
    return DiscreteSkyMap(skymap_config)


def initialize_butler(
    repo_path: Path,
    euclid_data_base_path: Path,
    tracts_from_dataset: str = "coadd",
    read_skymap_config: bool = False,
):
    ##===== Check that it's safe to start a new butler.
    repo_yaml = repo_path / "butler.yaml"
    if repo_yaml.exists():
        msg = f"Existing {repo_path}:\n    Cannot re-register over existing butler."
        raise FileExistsError(msg)

    ##===== Get the skymap (either read existing, or build new)
    if read_skymap_config:
        skymap_path = repo_path / euclid_constants.EUCLID_SKYMAP_FILENAME
        skymap = DiscreteSkyMap(skymap_path)
    else:
        skymap = build_euclid_skymap(
            repo_path, euclid_data_base_path, tracts_from_dataset=tracts_from_dataset
        )

    tract = skymap[0]

    logger.info(f"sky map has {tract.num_patches} patches")
    logger.info(f"example tract has bbox:\n    {tract._bbox}")

    idx = Index2D(tract.num_patches.x - 1, tract.num_patches.y - 1)
    patch = tract[idx]
    logger.info(f"example patch has inner_bbox:\n    {patch.inner_bbox}")
    # print(patch.patchInnerDimensions)

    ##===== Make the butler and register the skymap
    seed_config = {}
    if seed_config:
        config_str = json.dumps(seed_config, indent=2)
        logger.info(f"use seed config {config_str}")
    else:
        logger.info("use empty seed config")

    config = Config(seed_config)  # Must use Config, not ButlerConfig for make repo...?
    # config["datastore.formatters"] = ... # If custom formatters needed.
    butler_config = Butler.makeRepo(repo_path, config=config)  # static method

    # Watch below for 'writable' vs 'writeable' (extra 'e')
    butler = Butler.from_config(repo_path, writeable=True)  # classmethod preferred.
    skymap.register(euclid_constants.EUCLID_SKYMAP_NAME, butler)

    ###===== Register the euclid instrument - points to new 'ExposureL' StorageClass.
    EuclidInstrument().register(butler.registry)

    ##===== Register new DatasetTypes and their Formatters
    register_dataset_types(butler)  # Add the info for the coadds.

    ###===== Register the "main" run ('Euclid/Q1')
    butler.collections.register(euclid_constants.EUCLID_COLLECTION_NAME, CollectionType.CHAINED)


def extract_patch_data_from_fits_file(
    filepath: Path, skymap: DiscreteSkyMap, hdu: int, bandpass_regex: str = None
) -> list[tuple[str, dict]]:

    if bandpass_regex:
        # This is a bit gross - but the filter info for eg. the grid-psf images is not
        # in the file header. So we have to do a bit of regex to strip them out...
        # The regex pattern is slightly different for catalog-psf and grid-psf too, so
        # we switch them out.
        pattern = re.compile(bandpass_regex)
        matches = pattern.match(filepath.name)
        if not matches:
            msg = f"No matches in regex\n    {bandpass_regex}\n {filepath.stem}"
            raise ValueError(msg)
        regex_band = matches.group(1)  # Is 1-indexed for some reason...
    else:
        regex_band = None

    with fits.open(filepath) as f:
        header = f[hdu].header
        im_wcs = WCS(header)
        im_shape = [header.get("NAXIS1"), header.get("NAXIS2")]
        center_pix = (header.get("NAXIS1") // 2, header.get("NAXIS2") // 2)

        band = None
        if "FILTER" in header.keys():
            physical_filter = header["FILTER"]
            band = PHYSICAL_FILTER_TO_BAND[physical_filter]
        else:
            band = regex_band or None

    center_coord = im_wcs.pixel_to_world(*center_pix)
    center_point = SpherePoint(
        center_coord.ra.deg * degrees, center_coord.dec.deg * degrees
    )

    tract: TractInfo = skymap.findTract(center_point)
    # patch: PatchInfo = tract.findPatch(center_point)

    nx, ny = tract.getNumPatches()
    file_data = []
    for ii, patch in enumerate(tract):
        patch_data = dict(
            skymap=euclid_constants.EUCLID_SKYMAP_NAME,
            tract=tract.getId(),
            patch=patch.getSequentialIndex(),
            band=band,
            filepath=filepath,
        )
        file_data.append(patch_data)

    return file_data

def ingest_dataset(
    dataset: str, 
    butler: Butler, 
    euclid_data_base_path: Path, 
    cameras=("VIS", "NISP"), # Maybe in the future we want to use DECAM too?
):

    dataset_filelist = collect_files(dataset, euclid_data_base_path, cameras=cameras)
    hdu = HDU_LOOKUP[dataset]
    dataset_type = dataset_type_lookup[dataset]
    collection_name = collection_name_lookup[dataset]

    # Register dataset in a run HERE
    # butler "ingest" calls seem to auto-register collections, but "put" do not.
    logger.info(f"register new RUN collection: {collection_name}")
    butler.collections.register(
        collection_name, CollectionType.RUN
    )  # idempotent operation anyway.

    # Get the existing chain, append this new name, set the updated chain
    chained_run_name = euclid_constants.EUCLID_COLLECTION_NAME
    existing_chain = list(butler.registry.getCollectionChain(chained_run_name))
    logger.info(f"update existing run chain:\n    {existing_chain}")

    butler.registry.setCollectionChain(
        chained_run_name, existing_chain + [collection_name]
    )
    new_chain = list(butler.registry.getCollectionChain(chained_run_name))
    logger.info(f"new chain:\n    {new_chain}")

    skymap = butler.get(
        "skyMap",
        dataId={"skymap": euclid_constants.EUCLID_SKYMAP_NAME},
        collections="skymaps",
    )

    manifest: list[dict] = []

    # Need to extract filter from filename for dataset grid-psf.
    bandpass_regex = BANDPASS_REGEX_LOOKUP.get(dataset, None)

    converter = PUT_CONVERTER_LOOKUP.get(dataset, None)  
    if converter is not None:
        logger.info(f"use {converter.__name__} to convert for butler PUT")
    else:
        logger.info("use symlink with Butler INGEST")


    for ii, filepath in tqdm.tqdm(
        enumerate(dataset_filelist), 
        position=1, 
        desc="get dataId in FITSes", 
        total=len(dataset_filelist),
        leave=False
    ):
        file_data = extract_patch_data_from_fits_file(
            filepath, skymap, hdu, bandpass_regex=bandpass_regex
        )
        manifest.extend(file_data)

    manifest_table = Table(rows=manifest)
    MANIFEST_TABLES_PATH.mkdir(exist_ok=True)
    manifest_table_filepath = MANIFEST_TABLES_PATH / f"{dataset}_manifest.csv"
    manifest_table.write(manifest_table_filepath, overwrite=True)

    #manifest_table = manifest_table[manifest_table["tract"] == 202]

    dataset_type_obj = butler.registry.getDatasetType(dataset_type)
    tract_grouped = manifest_table.group_by("tract")
    for tract_group in tqdm.tqdm(tract_grouped.groups, position=1, desc="tracts"):
        tract_group.sort("patch")  # So that patch 0 is written first.

        if converter is not None:
            row0, remaining_rows = tract_group[0], tract_group[1:]

            dataset_ref = _make_dataset_ref(
                row0, dataset, dataset_type_obj, collection_name, butler
            )
            # Physically write the first patch to disk in `NumpyArray`.
            physical_uri = _write_physical_dataset(
                row0["filepath"], dataset_ref, butler, converter
            )
            # Symlink the rest.
            # for row in tqdm.tqdm(
            #     remaining_rows, position=2, desc="patches", leave=False
            # ):
            #     dataset_ref = _make_dataset_ref(
            #         row, dataset, dataset_type_obj, collection_name, butler
            #     )
            #     _symlink_dataset(physical_uri, dataset_ref, butler)

        else:
            for row in tqdm.tqdm(
                tract_group, position=2, desc="tract files", leave=False
            ):
                dataset_ref = _make_dataset_ref(
                    row, dataset, dataset_type_obj, collection_name, butler
                )
                _symlink_dataset(row["filepath"], dataset_ref, butler)


def _standardize_data_coord(row, dataset: str, butler: Butler):
    data_coord = dict(skymap=row["skymap"], tract=row["tract"], patch=row["patch"])
    if dataset != "segmap":
        data_coord["band"] = row["band"]
    return DataCoordinate.standardize(**data_coord, universe=butler.dimensions)


def _make_dataset_ref(
    row, dataset: str, dataset_type: str, collection_name: str, butler: Butler
):
    dataId = _standardize_data_coord(row, dataset, butler)
    return DatasetRef(datasetType=dataset_type, dataId=dataId, run=collection_name)


def _write_physical_dataset(
    filepath: Path, dataset_ref: DatasetRef, butler: Butler, converter: Callable
):
    """Convert the FITS segmap to a pickle once and put() it."""
    data = converter(filepath)
    physical_ref = butler.put(data, dataset_ref)
    return butler.getURI(physical_ref)


def _symlink_dataset(filepath: Path, dataset_ref: DatasetRef, butler: Butler):
    butler.ingest(
        FileDataset(path=filepath, refs=[dataset_ref]),
        transfer="symlink",
    )


def logging_setup(log_dir: Path = None, filename: str = "out.log"):

    handlers = [logging.StreamHandler()]  # []
    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(exist_ok=True, parents=True)
        log_filepath = log_dir / filename
        file_handler = TimedRotatingFileHandler(log_filepath, when="midnight")
        handlers.append(file_handler)

    logging.basicConfig(
        handlers=handlers,
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%y-%m-%d %H:%M:%S",
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--init", action="store_true", default=False)
    parser.add_argument(
        "-d", "--dataset", default=None, choices=ALLOWED_DATASETS, nargs="+"
    )
    parser.add_argument("-c", "--cameras", choices=ALLOWED_CAMERA_SUBDIRS)
    args = parser.parse_args()

    logging_setup()

    FILE_DATA_PATH = os.environ.get("EUCLID_FILE_DATA_PATH", None)
    # RAW_DATA_PATH = Path("/lustre/idac/plain/db/raw/euclid/q1")

    BUTLER_REPO = os.environ.get("EUCLID_BUTLER_REPO", None)
    if BUTLER_REPO is None or FILE_DATA_PATH is None:
        print(
            "load butler-config first:\n    "
            "it will set env variables EUCLID_FILE_DATA_PATH and EUCLID_BUTLER_REPO\n    "
            "\033[36;1msource /path/to/obs_euclid/butler-config\033[0m\nExiting."
        )
        sys.exit()

    FILE_DATA_PATH = Path(FILE_DATA_PATH)
    BUTLER_REPO = Path(BUTLER_REPO)
    if args.init:
        BUTLER_REPO.mkdir(exist_ok=True)
        initialize_butler(BUTLER_REPO, FILE_DATA_PATH)
        print(f"New bulter initialized in\n    \033[0;1m{BUTLER_REPO}\033[0m")
        sys.exit()

    butler_config = BUTLER_REPO / "butler.yaml"
    if not butler_config.exists():
        print("First initialize the butler \033[36;1mwith --init\033[0m")
        sys.exit()

    print(f"read Euclid data from:\n    \033[36;1m{FILE_DATA_PATH}\033[0m")

    if not args.dataset:
        msg = ", ".join(f"'{im}'" for im in ALLOWED_DATASETS)
        print(f"Provide dataset(s) with -d / --dataset from:\n    {msg}")
        sys.exit()

    # skymap = read_euclid_q1_skymap(BUTLER_REPO)
    butler = Butler(BUTLER_REPO, writeable=True)  # watch 'writeable' vs 'writable'

    dataset_bar_generator = tqdm.tqdm(args.dataset, position=0)
    for dataset in dataset_bar_generator:
        dataset_bar_generator.set_description(f"dataset='{dataset}'")
        ingest_dataset(dataset, butler, FILE_DATA_PATH)
