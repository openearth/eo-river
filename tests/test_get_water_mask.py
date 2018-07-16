import pytest

from eo_river import algorithms


def test_get_water_mask():
    region = {
        "geodesic": False,
        "type": "Polygon",
        "coordinates": [[
            [5.986862182617186, 52.517369933821186],
            [6.030635833740234, 52.517369933821186],
            [6.030635833740234, 52.535439735112924],
            [5.986862182617186, 52.535439735112924],
            [5.986862182617186, 52.517369933821186]
        ]]
    }

    start = '2017-01-01'
    stop = '2017-06-01'

    # these are also defaults
    percentile = 10
    ndwi_threshold = 0,
    scale = 10
    water_mask = algorithms.get_water_mask(region, start, stop,
                                           percentile, ndwi_threshold, scale)

    print(water_mask)

    # assert len(water_mask['features']) > 0
