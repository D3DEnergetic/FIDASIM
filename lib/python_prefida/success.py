#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.colored import colored


def success(string):
    #+##`success, str`
    #+Print a success message
    #+###Arguments
    #+     **str**: message
    #+
    #+###Example Usage
    #+```idl
    #+IDL> success, "Yay!!!"
    #+```
    print(colored('SUCCESS: ' + string, 'green'))
###############################################################################
if __name__ == "__main__":
    success('yay')