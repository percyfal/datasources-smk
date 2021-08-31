import os
import re
import yaml
import urllib
import numpy as np
import glob
from snakemake.utils import validate
import pandas as pd
from url import UrlMap

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
WRAPPER_PREFIX = workflow.wrapper_prefix
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX_RAW:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/percyfal/datasources-smk/main/workflow/wrappers"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

def _read(infile, index, schema, idcols=None):
    if infile is None:
        return None
    if os.path.splitext(infile)[1] == ".yaml":
        with open(infile) as fh:
            data = yaml.load(fh, yaml.Loader)
        assert(isinstance(data, list))
        df = pd.DataFrame(data)
    elif os.path.splitext(infile)[1] == ".tsv":
        df = pd.read_csv(infile, sep="\t")
    df.set_index(index, drop=False, inplace=True)
    df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schema)
    return df

datasources = _read(config["datasources"], ["data"],
                    "../schemas/datasources.schema.yaml")

##################################################
## Formatting functions and other utilities
##################################################
def wildcards_or(items):
    items = [str(x) for x in set(items)]
    return f'({"|".join(items)})'

##################################################
# Input collection functions
##################################################
url_map = UrlMap(datasources.source.to_list(), datasources.data.to_list())


def all_input(wildcards):
    return url_map.uri_keys
