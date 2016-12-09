#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os


def get_version(fidasim_dir):
    #+#get_version
    #+ Gets FIDASIM version number from git.
    #+ Falls back to reading VERSION file when git is not available
    #+***
    #+##Input Arguments
    #+    **fidasim_dir**: FIDASIM install directory
    #+
    #+##Example Usage
    #+```idl
    #+IDL> version = get_version(getenv("FIDASIM_DIR"))
    #+```
    version = ''

    ver_file = '{}{}VERSION'.format(fidasim_dir, os.path.sep)

    if os.path.isfile(ver_file):
        with open(ver_file) as f:
            version = f.read()

    return version
###############################################################################
if __name__ == "__main__":
    from lib.python_prefida.get_fidasim_dir import get_fidasim_dir
    a = get_version(get_fidasim_dir())

    print(a)
