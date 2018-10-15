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

    start = '2010-01-01'
    stop = '2016-01-01'
    scale = 10
    crs = 'EPSG:3857'

    water_mask = algorithms.get_water_mask(region, start, stop, scale, crs)

    print(water_mask)

    assert len(water_mask['features']) > 0


def test_get_water_network():
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

    start = '2010-01-01'
    stop = '2016-01-01'
    scale = 10
    crs = 'EPSG:3857'

    network = algorithms.get_water_network(region, start, stop, scale, crs)

    print(network)

    assert len(network['features']) > 0


def test_get_water_network_properties():
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

    start = '2010-01-01'
    stop = '2016-01-01'
    scale = 10
    step = 100
    crs = 'EPSG:3857'

    features = algorithms.get_water_network_properties(region, start, stop,
                                                       scale, crs, step)

    print(features)

    assert len(features['features']) > 0
