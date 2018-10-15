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
def get_water_mask(region, output, filter_bounds, start, stop):
    click.echo("Generating water mask from satellite data for %s - %s ..." %
               (start, stop))

    # read region GeoJSON file
    region = geojson.loads(open(region, 'r').read())

    # query water mask
    water_mask = eo_river.algorithms.get_water_mask(region, start, stop)

    # write results
    with open(output, 'w') as f:
        f.write(water_mask)


@click.command(name="get-network")
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
def get_network(region, output, filter_bounds, start, stop):
    click.echo("Generating water mask network from satellite data for %s - "
               "%s ..." % (start, stop))

    # read region GeoJSON file
    region = geojson.loads(open(region, 'r').read())

    # query water mask network
    network = eo_river.algorithms.get_water_network(region, start, stop)

    # write results
    with open(output, 'w') as f:
        f.write(network)


cli.add_command(get_water_mask)
cli.add_command(get_network)

if __name__ == "__main__":
    sys.exit(cli())  # pragma: no cover
