#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.colored import colored


def info(string):
    #+#info
    #+Print a informational message
    #+***
    #+##Arguments
    #+     **str**: message
    #+
    #+##Example Usage
    #+```python
    #+>>> info("This is an informative message")
    #+```
    print(colored('INFO: ' + string, 'cyan'))
