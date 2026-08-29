# obs_euclid

A set of simple scripts that work with the LSST stack to use
the an lsst.daf.butler.Butler with the Euclid MER data products.

You will need an installed copy of the lsst-stack.

## Prereqs

In your working dir `wkdir`, clone this repo:

`git clone git@github.com:aidansedgewick/obs_euclid.git`

You need an installed copy of the `lsst-scipipe` stack.

If you don't have one, install it. You don't need root permissions.
- Make a sensible new directory (eg. `<wkdir>/lsst_stack`). I don't think it's possible to move this later.
- `cd lsst_stack`
- Download the installer: `curl -OL https://ls.st/lsstinstall`
- Make the installer executable. `chmod 755 lsstinstall`
- Install and specify a version: `./lsstinstall -T v30_0_10`.
  - This will probably take 5/10 mins.
  - You can find what versions are new here: https://rsp.lsst.io/updates/index.html
- `source loadLSST.sh`
- `eups distrib install -t v30_0_10 lsst_distrib`  
  - This step also probably takes 10 mins.

Source & setup your new `lsst-scipipe`

You can then install the new `obs_euclid` on top:
- `cd <wkdir>/obs_euclid`
- `python3 -m pip intstall -e .`
  - Install in editable mode, in case you need to change something.

## Data layout.

The FITS file data are organised as they are downloaded from the IRSA servers
https://irsa.ipac.caltech.edu/ibe/data/euclid/q1/

ie.
```
/path/to/data/euclid/data/
  MER/
    [tile-id]/
      NISP/
        EUC_MER_{datasets}_NIR-{H,J,Y}_[tile-id]_[rng]_[isot].fits
      VIS/
        EUC_MER_{datasets}_VIS_[tile-id]_[rng]_[isot].fits
  MER_SEG/
    [tile-id]/
      EUC_MER_FINAL_SEGMAP_[].fits
```

where datasets are:
  - BGMOD (ingest with key `bgmod`)
  - BGSUB-MOSAIC (ingest with key `coadd`)
  - CATALOG-PSF (`catalog-psf`)
  - GRID-PSF (`grid-psf`)
  - MOSAIC-*-RMS (`rms`)
  - MOSAIC-*-FLAG (`flag`)

and separately:
  - SEGMAP (`segmap`)

## Creating the butler

### Setup

To create the butlerized data:

Add two variables to a file `q1-butler-setup`
```
export EUCLID_FILE_DATA_PATH=/path/to/euclid/data
export EUCLID_BUTLER_REPO=/path/to/repo
```
and source it.

The names for the skymap and dataset types and collections are defined in `lsst.obs.euclid.euclid_constants`

### Initializing the butler.

Initialize the new butler with:
  `python3 obs_euclid/scripts/auto_butler.py --init`

It does two main things: creates a `skyMap`, and registers all of the `datasetTypes`

#### Skymap

The skymap is called `euclid_q1`.

Initializing will scrape through all of the VIS coadds and extract the WCS information
and build an `lsst.skymap.DiscreteSkyMap`, where each tract is defined as a single
Euclid tile, as one very large patch.

There are 24 of 352 tiles that are not the "normal" 19200x19200 pixels.
The largest are (21600, 19200) and (19200, 21600) - so max 21600 in both x/y
So actually we set the patch size in pixels to (21600, 21600).
The Butler is flexible enough that when it loads the data from FITS file,
it will read the all of the data.

(If we set the patchsize = (19200, 19200) pixels, I think it would force
the oversized tiles to be actually four patches one "normal size" [0:19200,0:19200]
two long and thin [0:19199,19200:x_size] and [19200:y_size,0:19199] and one very small
[19200:y_size,19200:x_size]. Or the three "extra" patches would simply be sliced off.)

#### Dateset types

The dataset types are defined/registered after skymap.
They are defined & registered in a module `lsst.obs.euclid.define_datasets`

##

There are a few datasets/collections which are the `ExposureF` (float32) storage class.

- `euclidCoadd`
- `euclidBgMod`
- `euclidGridPsf`
- `euclidRmsMap`

They are all ingested as `symlink`, so no actual data is written for each of these files.
They are loaded in about 2 mins each.

A quirk with the `grid-psf` (`euclidGridPsf` collection): as the data for these are stored in HDU1, not HDU0, there is some error on load.
But the actual data seems like it loads correctly.

##

There are also
- `euclidFlagMap`
are storage class `ExposureI` (int32).

They also seem to load and read fine.

##

I have some issues with the `euclidCatalogPsf` catalog-psf datasets.
They must be physically written to disk in a format that the Butler can then re-read.
I think that I am reading the WCS from the FITS headers incorrectly, because several catalog-psf files are mapped to the same tract.

##

The segmentation maps (`euclidSegMap`) are also a bit problematic. 
Because these files are int64 FITS, I think there is no native Butler storage class that can them.

The smart solution by from Jennifer Li at UIUC is to read these in and write them to the `NumpyArray` storage class.
However this means that they need to be physically written to disk.
This took around 30mins using the hardware at DK-IDAC.
It also means that there is unfornately no WCS information stored with them. 

But as the WCS for all image products is the same for a single tile, when you find the
bounding box for a product with WCS available (eg. coadds), the same bbox should apply to
the segmap...


### Loading the data

Load the datasets with eg.:
-  `python3 obs_euclid/scripts/auto_butler.py --dataset coadd`
-  `python3 obs_euclid/scripts/auto_butler.py --dataset rms`
-  `python3 obs_euclid/scripts/auto_butler.py --dataset bgmod`
-  `python3 obs_euclid/scripts/auto_butler.py --dataset flag`

The Chained collection is called `Euclid/Q1`.

You can see/change all of these names in 

Some of the products ingest very well.

eg. coadd, rms, bgmod, flag.


## Querying the butler

You can make basic queries with Butler.query_datasets().
```
>>> from lsst.daf.butler import Butler
>>> b = Butler("/path/to/repo")
>>> coadd_refs = b.query_datasets(
...    "euclidCoadd", where="band.name='J', limit=10, collections="Euclid/Q1"
... )
>>> print(coadd_refs[0])
>>> print(coadd_refs[0].ospath)
>>> im_data = b.get(coadd_refs[0])
>>> print(im.wcs) # wcs to make transforms
>>> print(im.image.array) # pixel data
```

Because the skymap is made of tracts, you can do spatial queries:
```
>>> ra = 60.3132
>>> dec = -51.0
>>> spatial_q = f"band.name='J' AND tract.region OVERLAPS POINT({ra}, {dec})"
>>> spatial_refs = b.query_dataset(
...     "euclidCoadd", where=spatial_q, collections="Euclid/Q1"
... )
>>> print(len(spatial_refs)) # Should be one, at the centre of a tile
>>> im = b.get(spatial_refs[0])
```

## Issues

- non-standard tract size means that tract WCS does not exactly match with loaded image WCS.
  - currently you have to take cutouts from the full loaded image, 
    rather than requesting a windowed read (that is, you can't currently provide 'bbox to `Butler.get()`).
- catalog-psf tables are not ingested correctly.
- Ingesting images segmaps without physically writing to disk - would give faster ingest.
