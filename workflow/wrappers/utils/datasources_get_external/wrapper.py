#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

source = snakemake.input.get("uri", None)
if source is None or len(source) == 0:
    source = snakemake.params.uri
target = snakemake.output.target
scheme = snakemake.params.scheme

cmd_map = {"": "ln -s", "rsync": "rsync -av", "sftp": "cp", "file":
           "ln -s", "http": "wget", "https": "wget"}
cmd = cmd_map[scheme]

outdir = os.path.dirname(target)
shell("mkdir -p {outdir}")

if cmd == "ln -s":
    source = os.path.abspath(source)

logger.debug(f"Using command {cmd} {source} {target} {log}")

if cmd == "wget":
    shell("{cmd} {source} -O {target} {log}")
else:
    shell("{cmd} {source} {target} {log}")
