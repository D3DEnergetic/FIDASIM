#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.colored import colored


def warn(string):
    #+##`warn, string`
    #+Print a warning message
    #+###Arguments
    #+     **string**: message
    #+
    #+###Example Usage
    #+```idl
    #+IDL> warn, "This may be a problem"
    #+```
    print(colored('WARNING: ' + string, 'magenta'))
###############################################################################
if __name__ == "__main__":
    warn('test')