#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    print(colored('ERROR: {}'.format(string), c='r'))

    if halt:
        raise Exception(string)
