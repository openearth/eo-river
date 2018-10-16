# -*- coding: utf-8 -*-

"""Console script for eo-river."""
import sys
import click

import geojson

import eo_river.algorithms


@click.group()
def cli():
    pass


@click.command(name="get-water-mask")
@click.option("-r", "--region", required=True,
              help="Path of the input region as GeoJSON file.")
@click.option("-o", "--output", required=True,
              help="Output water mask file path as GeoJSON.")
@click.option("-f", "--filter-bounds", required=False,
              help="Path of the input GeoJSON file defining intersecting "
                   "geometries to be filtered.")
@click.option("--start", default="2010-01-01",
              help="Start time (default is 2010-01-01).")
@click.option("--stop", default="2015-01-01",
              help="Stop time (default is 2015-01-01).")
@click.option("--scale", default=10,
              help="Scale of the figure (default is 10)")
@click.option("--crs", default="EPSG:4326",
              help="Coordinate system as an EPSG code (default is 'EPSG:4326').")
def get_water_mask(region, output, filter_bounds, start, stop, scale, crs):
    click.echo("Generating water mask from satellite data for %s - %s ..." %
               (start, stop))

    # read region GeoJSON file
    region = geojson.loads(open(region, 'r').read())

    # query water mask
    water_mask = eo_river.algorithms.get_water_mask(region, start, stop, scale, crs)

    # write results
    with open(output, 'w') as f:
        f.write(geojson.dumps(water_mask))


@click.command(name="get-water-network")
@click.option("-r", "--region", required=True,
              help="Path of the input region as GeoJSON file.")
@click.option("-o", "--output", required=True,
              help="Output netwerk file path as GeoJSON.")
@click.option("-f", "--filter-bounds", required=False,
              help="Path of the input GeoJSON file defining intersecting "
                   "geometries to be filtered.")
@click.option("--start", default="2010-01-01",
              help="Start time (default is 2010-01-01).")
@click.option("--stop", default="2015-01-01",
              help="Stop time (default is 2015-01-01).")
@click.option("--scale", default=10,
              help="Scale of the figure (default is 10)")
@click.option("--crs", default="EPSG:4326",
              help="Coordinate system as an EPSG code (default is 'EPSG:4326').")
def get_network(region, output, filter_bounds, start, stop, scale, crs):
    click.echo("Generating water mask network from satellite data for %s - "
               "%s ..." % (start, stop))

    # read region GeoJSON file
    region = geojson.loads(open(region, 'r').read())

    # query water mask network
    network = eo_river.algorithms.get_water_network(region, start, stop, scale, crs)

    # write results
    with open(output, 'w') as f:
        f.write(geojson.dumps(network))

        
@click.command(name="get-water-network-properties")
@click.option("-r", "--region", required=True,
              help="Path of the input region as GeoJSON file.")
@click.option("-o", "--output", required=True,
              help="Output network file path as GeoJSON.")
@click.option("-f", "--filter-bounds", required=False,
              help="Path of the input GeoJSON file defining intersecting "
                   "geometries to be filtered.")
@click.option("--start", default="2010-01-01",
              help="Start time (default is 2010-01-01).")
@click.option("--stop", default="2015-01-01",
              help="Stop time (default is 2015-01-01).")
@click.option("--scale", default=10,
              help="Scale of the figure (default is 10)")
@click.option("--crs", default="EPSG:4326",
              help="Coordinate system as an EPSG code (default is 'EPSG:4326').")
@click.option("--step", default=100,
              help="Distance between features (default is 100).")
def get_network_properties(region, output, filter_bounds, start, stop, scale, crs, step):
    click.echo("Generating water mask network including properies from satellite data for %s - "
               "%s ..." % (start, stop))

    # read region GeoJSON file
    region = geojson.loads(open(region, 'r').read())

    # query water mask network with properties
    features = eo_river.algorithms.get_water_network_properties(region, start, stop, scale, crs, step)

    # write results
    with open(output, 'w') as f:
        f.write(geojson.dumps(features))


cli.add_command(get_water_mask)
cli.add_command(get_network)
cli.add_command(get_network_properties)

if __name__ == "__main__":
    sys.exit(cli())  # pragma: no cover
