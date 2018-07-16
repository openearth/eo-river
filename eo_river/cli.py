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


# @click.command(name="get-water-mask")
# @click.option("-o", "--options-file", required=True, help="Options file in YAML format")
# @click.option("-r", "--results-dir", required=True, help="Result directory")
# @click.options("-l", "--local", required=False, help="Use locally installed generators")
# def generate_model(options_file, results_dir):
#     # two YAML docs are expected in this file, one generic and one model specific
#     dicts = builder.parse_config(options_file)
#     genopt, modopt = dicts
#     # TODO validate config
#     # TODO fill in all defaults (for now we should supply all)
#     msg = f"Going to create a '{genopt['model']}'/'{modopt['concept']}' model, it will be placed in '{results_dir}'"
#     general_options(genopt)
#     click.echo(msg)
#     click.echo(builder.get_generator_names())



if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
