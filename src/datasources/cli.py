"""Console script for datasources"""
import logging
import sys
from argparse import ArgumentParser

from . import __version__

__author__ = "Per Unneberg"


logger = logging.getLogger(__name__)


def main(arg_list=None):
    if arg_list is None:
        arg_list = sys.argv[1:]
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = ArgumentParser(description=__doc__, prog="datasources")
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    parser.add_argument(
        "--debug", action="store_true", default=False, help="Print debug messages"
    )

    args = parser.parse_args(arg_list)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    return 0
