import os
import re
import sys
import yaml
import urllib
import numpy as np
import glob
from snakemake.utils import validate
import pandas as pd

try:
    from datasources.url import UrlMap
    from datasources.utils import wildcards_or
except ModuleNotFoundError:
    print("datasources module not installed")
    print("Either: ")
    print("  1. point PYTHONPATH to workflow source directory or")
    print("  2. run 'python -m pip install -e /path/to/datasources' or")
    print(
        "  3. run 'python -m pip install git+https://github.com/percyfal/datasources-smk@main'"
    )
    raise


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


# ##### load datasources #####


def _read(infile, index, schema, idcols=None):
    if infile is None:
        return None
    if os.path.splitext(infile)[1] == ".yaml":
        with open(infile) as fh:
            data = yaml.load(fh, yaml.Loader)
        assert isinstance(data, list)
        df = pd.DataFrame(data)
    elif os.path.splitext(infile)[1] == ".tsv":
        df = pd.read_csv(infile, sep="\t")
    df.set_index(index, drop=False, inplace=True)
    df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schema)
    return df


if config.get("datasources") is None:
    configfiles = [
        "datasources.yaml",
        "datasources.tsv",
        "config/datasources.yaml",
        "config/datasources.tsv",
    ]
else:
    configfiles = [config.get("datasources")]

datasources = None

for infile in configfiles:
    if infile is not None and os.path.exists(infile):
        datasources = _read(infile, ["data"], "../schemas/datasources.schema.yaml")
    if datasources is not None:
        break
if datasources is None:
    print("No datasources found: exiting")
    sys.exit(1)


##################################################
# Input collection functions
##################################################
url_map = UrlMap(datasources.source.to_list(), datasources.data.to_list())


def all_input(wildcards):
    return url_map.uri_keys
