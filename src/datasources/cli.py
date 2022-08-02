"""Console script for datasources"""
import logging
import sys
import subprocess as sp
from argparse import ArgumentParser

from . import __version__

__author__ = "Per Unneberg"


logger = logging.getLogger(__name__)


def run(args):
    print(f"snakemake {' '.join(args.extra_options)}")



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
        help=(
            "Run snakemake workflow"
        ),
    )
    parser.set_defaults(runner=run)

    args, extra = top_parser.parse_known_args(arg_list)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    args.extra_options = extra

    args.runner(args)

    return 0
