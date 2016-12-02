#!/usr/bin/env python
# -*- coding: utf-8 -*-

from python_prefida import colored


def warn(string):
    """
    #+##`warn, string`
    #+Print a warning message
    #+###Arguments
    #+     **string**: message
    #+
    #+###Example Usage
    #+```idl
    #+IDL> warn, "This may be a problem"
    #+```
    """
    print(colored('WARNING: ' + string, c='y'))
