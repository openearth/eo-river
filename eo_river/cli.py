# -*- coding: utf-8 -*-

"""Console script for eo_river."""
import sys
import click

import geojson

import eo_river.algorithms


@click.group()
def main():
    return 0


@click.command(name="get-water-mask")
@click.option("-r", "--region", required=True,
              help="Path of the input region as GeoJSON file.")
@click.option("-r", "--output", required=True,
              help="Output water mask file path as GeoJSON.")
@click.option("--start", default="2017-01-01",
              help="Start time (default is 2017-01-01).")
@click.option("--stop", default="2017-06-01",
              help="Stop time (default is 2017-06-01).")
def get_water_mask(region, output, start, stop):
    click.echo("Generating water mask from satellite data for %s - %s ..." %
               (start, stop))

    # read region GeoJSON file
    region = geojson.loads(open(region, 'r').read())

    # query water mask
    water_mask = eo_river.algorithms.get_water_mask(region, start, stop)

    # write results
    f = open(output, 'w')
    f.write(water_mask)


main.add_command(get_water_mask)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
