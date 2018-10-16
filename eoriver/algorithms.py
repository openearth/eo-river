# -*- coding: utf-8 -*-

"""Main module."""

import hydroengine as he


def get_water_mask(region, start, stop, scale, crs):
    water_mask = he.get_water_mask(region, start, stop, scale, crs)

    return water_mask


def get_water_network(region, start, stop, scale, crs):
    water_network = he.get_water_network(region, start, stop, scale, crs)

    return water_network


def get_water_network_properties(region, start, stop, scale, step, crs):
    water_network_properties = he.get_water_network_properties(region, start,
                                                               stop, scale,
                                                               step, crs)

    return water_network_properties


def generate_network(watermask):
    """    This function loads a watermask, skeletonise and return as a network
    of branches.

    Args:
       watermask: geojson-string or -file with only the polygon of the river
       output: write output to file instead of returning output as string

    Output:
       network: geojson-string or -file with the network
    """
    import ratin.network as nw

    STARTEND = [[197727.048, 503291.164], [195687.574,
                                           504118.171]]  # This should be derived automatically!

    SMOOTHDEGREE, cutpoints, densi, space, horiLines = 15.0, 20, 5.0, 5.0, 15  # dense grid

    # This now works for single branch
    network = Network()
    network.load_geometries(polygons=watermask, se_points=STARTEND)
    network.densify(densi)
    network.construct(spacing=[space])
    network.delete(branch='1', num_of_pnts=cutpoints)
    # Representation
    network.plot(1)
    # Network Stats:
    network.overview(branch='1', fignum=2, smoothdegree=SMOOTHDEGREE,
                     printoutput=False, units='radians')
    # Save output of network file
    # network.save('output/network')

    return None
