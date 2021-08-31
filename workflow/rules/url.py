#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2021, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import re
from os.path import dirname, basename, join as pjoin
import glob
from snakemake import WorkflowError, logger
from urllib import parse


# Constructor-like wrapper
def urlparse(url, scheme=None):
    if url is None:
        return None
    if not scheme is None:
        if isinstance(url, str):
            url = f"{scheme}:{url}"
        elif isinstance(url, list):
            url = [f"{scheme}:{u}" for u in url]
    if isinstance(url, list):
        return [Url(**parse.urlparse(u)._asdict()) for u in url]
    return Url(**parse.urlparse(url)._asdict())


# Don't call directly; use wrapper
class Url(parse.ParseResult):
    wc_regex = re.compile(r"[*\]\]?]+")

    def __init__(self, *args, **kwargs):
        super().__init__()


    def _uri(self, scheme=False):
        if scheme:
            scheme = self.scheme
        else:
            scheme = ''
        uri = parse.urlunsplit((scheme, self.netloc, self.path, '', ''))
        uri = re.sub("^//", "", uri)
        return uri


    @property
    def uri(self):
        """Format url as uri without scheme"""
        return self._uri()


    @property
    def full_uri(self):
        """Format url as uri with scheme"""
        return self._uri(scheme=True)


    @property
    def source(self):
        """Create snakemake target representation of data source"""
        uri = self.uri
        if self.scheme == 'sftp':
            # FIXME: this setup currently fails; in contrast to fabric
            # it seems paramiko doesn't make use of the ssh-agent
            try:
                from snakemake.remote.SFTP import RemoteProvider
                SFTP = RemoteProvider()
                uri = SFTP.remote(self.uri)
            except WorkflowError as e:
                logger.error(e)
        elif self.scheme == 'http' or self.scheme == 'https':
            try:
                from snakemake.remote.HTTP import RemoteProvider
                HTTP = RemoteProvider()
                uri = HTTP.remote(self.uri)
            except WorkflowError as e:
                logger.error(e)
        return uri


    @property
    def wildcard(self):
        if self.wc_regex.search(self.uri) is None:
            return None
        return self.basename


    def glob(self):
        if self.netloc == '':
            return urlparse(glob.glob(self.uri), scheme=self.scheme)
        try:
            import fabric
            import fnmatch
        except ImportError:
            print("modules 'fabric' and 'fnmatch' needed for remote wildcards")
            raise
        try:
            c = fabric.Connection(self.netloc)
        except Exception as e:
            print(e)
            raise
        data = fnmatch.filter(c.sftp().listdir(dirname(self.path)), self.wildcard)
        ret = []
        for d in data:
            ret.append(parse.urlunsplit((self.scheme, self.netloc, pjoin(dirname(self.path), d), None, None)))
        return urlparse(ret)


    @property
    def basename(self):
        return basename(self.path)


    @property
    def dirname(self):
        return dirname(self.uri)


    def __str__(self):
        return self.uri




class UrlMap:
    def __init__(self, source, target):
        assert isinstance(source, list) and isinstance(target, list), \
            print("source and target must be list of strings")
        self._source = list()
        self._target = list()
        for (src, tgt) in zip(source, target):
            wc = urlparse(tgt).wildcard
            if wc is not None:
                tgt = urlparse(tgt)
                src = urlparse(pjoin(src, tgt.wildcard))
                globsrc = src.glob()
                self._source.extend(globsrc)
                self._target.extend(urlparse([pjoin(tgt.dirname, fn.basename) for fn in globsrc]))
            else:
                self._target.append(urlparse(tgt))
                self._source.append(urlparse(src))


    @property
    def source(self):
        return self._source


    @property
    def target(self):
        return self._target


    @property
    def uri_keys(self):
        return [x.uri for x in self.target]


    @property
    def uri_dict(self):
        return dict(zip(self.uri_keys, self.source))


    def get_source_url(self, wildcards):
        return self.uri_dict[wildcards.target]


    def get_source(self, wildcards):
        return self.get_source_url(wildcards).source


    def get_source_uri(self, wildcards):
        return self.get_source_url(wildcards).uri


    def get_source_scheme(self, wildcards):
        return self.get_source_url(wildcards).scheme


    @property
    def pairs(self):
        return zip(self.source, self.target)


    def __iter__(self):
        for (source, target) in self.pairs:
            yield source, target
