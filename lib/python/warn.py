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
    #+```python
    #+>>> warn("This may be a problem")
    #+```
    print(colored('WARNING: ' + string, 'magenta'))
