import pytest

from eo_river import eo_river

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

    water_mask = eo_river.get_water_mask(region)

    print(water_mask)

    assert len(water_mask['features']) > 0
