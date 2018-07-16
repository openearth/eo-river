from eo_river import eo_river


def test_generate_network():
    watermaskfile = r'test_resources\ijssel_N_polygon.shp'
    network = eo_river.generate_network(watermask=watermaskfile)

    assert(network)
