import os
import urllib
import numpy as np
from snakemake.utils import validate
import pandas as pd

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
WRAPPER_PREFIX = workflow.wrapper_prefix
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX_RAW:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/percifal/datasources-smk/main/workflow/wrappers"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

def _read_tsv(infile, index, schemafile):
    if infile is None:
        schema = _load_configfile(os.path.join(workflow.current_basedir, schemafile))
        cols = schema["required"]
        df = pd.DataFrame({k:[] for k in cols})
        df = df.set_index(index, drop=False)
    else:
        df = pd.read_csv(infile, sep="\t").set_index(index, drop=False)
        df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schemafile)
    ## Check for duplicate rows
    if any(df.duplicated()):
        logger.error(f"duplicated rows in {infile}:")
        logger.error(f"  {df[df.duplicated()]}")
        sys.exit(1)
    return df

datasources = _read_tsv(config["datasources"], ["data"],
                        "../schemas/datasources.schema.yaml")



##################################################
## Uri parsing functions
##################################################
def get_uri_scheme(uri):
    return urllib.parse.urlparse(uri).scheme

def get_uri_netloc(uri):
    return urllib.parse.urlparse(uri).netloc

def parse_uri(uri):
    """Parse uri and return snakemake target"""
    allowed_schemes = ['', 'rsync', 'file', 'sftp', 'http', 'https']
    scheme = get_uri_scheme(uri)
    uri = re.sub(f"{scheme}://", "", uri)
    if not scheme in allowed_schemes:
        logger.error(f"scheme '{scheme}' not allowed: use one of {','.join(allowed_schemes[1:])}")
        sys.exit(1)
    if scheme in ['', 'file'] and not uri.startswith("/"):
        uri = os.path.normpath(os.path.abspath(uri))
    if scheme == 'sftp':
        try:
            from snakemake.remote.SFTP import RemoteProvider
            SFTP = RemoteProvider()
            uri = SFTP.remote(uri)
        except WorkflowError as e:
            logger.error(e)
    if scheme == 'http' or scheme == 'https':
        try:
            from snakemake.remote.HTTP import RemoteProvider
            HTTP = RemoteProvider()
            uri = HTTP.remote(uri)
        except WorkflowError as e:
            logger.error(e)
    return uri

def datasources_get_external_input(uri):
    if get_uri_scheme(uri) == "sftp":
        return parse_uri(uri)
    netloc = get_uri_netloc(uri)
    uri = parse_uri(uri)
    if re.search("[:@]", netloc):
        return []
    return uri

##################################################
# Input collection functions
##################################################
def all_input(wildcards):
    d = datasources["data"].to_list()
    return d
