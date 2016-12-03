#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from lib.python_prefida import source_file


def get_fidasim_dir():
    """
    ;+#get_fidasim_dir
    ;+ Gets FIDASIM install directory
    ;+***
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> fida_dir = get_fidasim_dir()
    ;+```
    """
#    return file_dirname(file_dirname(source_file()))

    filepath = source_file()
    directory = os.path.dirname(os.path.dirname(filepath))
    return directory

###############################################################################
if __name__ == "__main__":
    print(get_fidasim_dir())
