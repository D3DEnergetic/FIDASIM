#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
from lib.python_prefida.source_file import source_file


def get_fidasim_dir():
    #+#get_fidasim_dir
    #+ Gets FIDASIM install directory
    #+***
    #+##Output Arguments
    #+     **directory**: FIDASIM install directory.
    #+##Example Usage
    #+```python
    #+>>> fida_dir = get_fidasim_dir()
    #+```
    filepath = source_file()

    directory = os.path.dirname(os.path.dirname(os.path.dirname(filepath)))

    return directory