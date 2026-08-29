from __future__ import annotations

from pathlib import Path

from lsst.daf.butler import StorageClassFactory
from lsst.obs.base import Instrument
from lsst.pipe.base import Pipeline
from lsst.utils import getPackageDir

# from lsst.utils.introspection import

from lsst.obs.euclid.euclid_filters import EUCLID_FILTER_DEFINITIONS

_PACKAGE_DIR = Path(__file__).parents[4]


class EuclidInstrument(Instrument):

    filterDefinitions = EUCLID_FILTER_DEFINITIONS

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def getName(cls) -> str:
        return "Euclid"

    @classmethod
    def getPackageDir(cls):
        return _PACKAGE_DIR

    def getCamera(self):
        raise NotImplementedError("No camera defined yet")

    def getRawFormatter(self):
        raise NotImplementedError("No raw formatter defined yet")

    def register(self, registry, update=False):

        instrument_records = {"name": self.getName()}
        with registry.transaction():
            registry.syncDimensionData(
                "instrument",
                instrument_records,
                update=update,
            )
            self._registerFilters(registry)

    # def _register_new_storage_classes(self):
    #     print("about to register int64 ExposureL")
    #     exposure_l_config_path = _PACKAGE_DIR / "config/storageClasses.yaml"

    #     sf = StorageClassFactory()
    #     sf.addFromConfig(exposure_l_config_path)

    @classmethod
    def applyConfigOverrides(cls, name, config):
        super().applyConfigOverrides(name, config)
        print("applying OVERRIDES NOW!")
        config["storageClasses.config"].append(
            str(_PACKAGE_DIR / "config/storageClasses.yaml")
        )
