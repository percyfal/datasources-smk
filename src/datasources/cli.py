"""Console script for datasources"""
import logging
import os
import subprocess as sp
import sys
from argparse import ArgumentParser

from . import __version__

__author__ = "Per Unneberg"


logger = logging.getLogger(__name__)

SNAKEFILE = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "workflow", "Snakefile"
)
__RUNDOC__ = (
    "Run datasources snakemake workflow"
    " "
    "Given a datasources configuration file, generate data files sources."
)


def run(args):
    datasources = ""
    if args.configfile is not None:
        datasources = f"--config datasources={args.configfile}"
    sp.check_output(
        f"snakemake -s {SNAKEFILE} {' '.join(args.extra_options)} {datasources}",
        shell=True,
    )


def main(arg_list=None):
    if arg_list is None:
        arg_list = sys.argv[1:]
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    top_parser = ArgumentParser(description=__doc__, prog="datasources")
    top_parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    top_parser.add_argument(
        "--debug", action="store_true", default=False, help="Print debug messages"
    )

    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    parser = subparsers.add_parser(
        "run",
        help="Run datasources snakemake workflow",
        description=__RUNDOC__,
    )
    parser.add_argument(
        "--configfile",
        action="store",
        default=None,
        help="datasources configuration file",
    )
    parser.set_defaults(runner=run)

    args, extra = top_parser.parse_known_args(arg_list)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    args.extra_options = extra

    args.runner(args)

    return 0
