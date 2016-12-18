#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.colored import colored


def error(string, halt=False):
    #+#error
    #+Print a error message
    #+***
    #+##Arguments
    #+     **string**: message
    #+
    #+##Keyword Arguments
    #+     **halt**: Halt program execution
    #+
    #+##Example Usage
    #+```python
    #+>>> error("Error message")
    #+```
    print(colored('ERROR: {}'.format(string), 'red'))

    if halt:
        raise Exception()
