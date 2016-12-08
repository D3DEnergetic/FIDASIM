#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib.python_prefida.colored import colored


def info(string):
    #+#info
    #+Print a informational message
    #+***
    #+##Arguments
    #+     **str**: message
    #+
    #+##Example Usage
    #+```idl
    #+IDL> info, "This is an informative message"
    #+```
    print(colored('INFO: ' + string, 'cyan'))
###############################################################################
if __name__ == "__main__":
    info('hello')