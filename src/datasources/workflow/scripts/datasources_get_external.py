#!/usr/bin/env python3
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

source = snakemake.input.get("uri", None)
if source is None or len(source) == 0:
    source = snakemake.params.uri
target = snakemake.output.target
scheme = snakemake.params.scheme

cmd_map = {
    "": "ln -s",
    "rsync": "rsync -av",
    "sftp": "cp",
    "file": "ln -s",
    "http": "wget",
    "https": "wget",
}
cmd = cmd_map[scheme]

if not isinstance(source, list):
    source = [source]
if not isinstance(target, list):
    target = [target]
if len(source) != len(target):
    print("source and target must be of equal lengths")
    raise Exception


# FIXME: Move to datasources function. Add pattern/tag argument for
# restricting outputs
for src, tgt in zip(source, target):
    outdir = os.path.dirname(tgt)
    shell("mkdir -p {outdir}")

    if cmd == "ln -s":
        src = os.path.abspath(src)

    logger.debug(f"Using command {cmd} {src} {tgt} {log}")

    if cmd == "wget":
        shell("{cmd} {src} -O {tgt} {log}")
    else:
        shell("{cmd} {src} {tgt} {log}")
