from lsst.obs.base import FilterDefinition, FilterDefinitionCollection

EUCLID_FILTER_DEFINITIONS = FilterDefinitionCollection(
    FilterDefinition(physical_filter="VIS", band="VIS"),
    FilterDefinition(physical_filter="NIR-Y", band="Y"),
    FilterDefinition(physical_filter="NIR-J", band="J"),
    FilterDefinition(physical_filter="NIR-H", band="H"),
)
