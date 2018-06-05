# -*- coding: utf-8 -*-

"""Console script for eo_river."""
import sys
import click


@click.command()
def main(args=None):
    """Console script for eo_river."""
    click.echo("Replace this message by putting your code into "
               "eo_river.cli.main")
    click.echo("See click documentation at http://click.pocoo.org/")
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
