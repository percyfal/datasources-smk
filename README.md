# Snakemake workflow: datasources-smk

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build status](https://github.com/percyfal/datasources-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/percyfal/datasources-smk/actions?query=workflow%3ATests) ![License](https://img.shields.io/badge/license-MIT-blue.svg)

## About

Snakemake workflow to setup external data for data analyses. The data
sources can be local or remote files.

## Installation

Easiest is to install via pip:

	python -m pip install git+https://github.com/percyfal/datasources-smk@main
	
Alternatively grab a copy of the source distribution and make a local
install:

	git clone https://github.com/percyfal/datasources-smk.git
	cd datasources-smk
	python -m pip install -e .

## Usage 

The workflow and additional commands run via the main entry point:

	datasources -h
	datasources run -j 1
	datasources run --configfile datasources.yaml

See the subcommand help for more information.


## Information

This workflow reads a datasources yaml file with list elements
consisting of `data` and `source` keys, or alternatively a
tab-separated file with columns `data` and `source`. The `data` and
`source` keys define file URI mapping from source to a snakemake
target. Supported URI schemes are currently `rsync`, `file`, `sftp`,
`http` and `https`. 

There are two optional keys; `description` is a free text field for
provenance information, and `tag` a tag to group data types such that
subsets of datasources can be targeted.

The datasources file can be provided via the `--configfile` option. If
unset, the workflow will look for files `datasources.yaml`,
`datasources.tsv`, `config/datasources.yaml` and
`config/datasources.tsv`, in that order.

URIs are given according to the [URI generic
syntax](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier#Syntax).
For instance, a local file is given as `file:relative/path/to/source`,
whereas examples of a remote files are
`rsync://example.com:80/absolute/path/to/source` and
`sftp://example.com:80/absolute/path/to/source`.


## Authors

* Per Unneberg (@percyfal)

## Testing

Test cases are in the subfolder `src/datasources/.test`. They are automatically
executed via continuous integration with [Github
Actions](https://github.com/features/actions).

