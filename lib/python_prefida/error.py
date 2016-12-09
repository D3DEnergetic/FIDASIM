#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.colored import colored


def error(string, halt=False):
    #+#error
    #+Print a error message
    #+***
    #+##Arguments
    #+     **str**: message
    #+
    #+##Keyword Arguments
    #+     **halt**: Halt program execution
    #+
    #+##Example Usage
    #+```idl
    #+IDL> error, "=("
    #+```
    print(colored('ERROR: {}'.format(string), 'red'))

    if halt:
        raise Exception()
###############################################################################
if __name__ == "__main__":
    error('wrong', halt=True)