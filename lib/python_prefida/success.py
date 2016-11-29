#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib import colored


def success(string):
    """
    #+##`success, str`
    #+Print a success message
    #+###Arguments
    #+     **str**: message
    #+
    #+###Example Usage
    #+```idl
    #+IDL> success, "Yay!!!"
    #+```
    """

    print colored('SUCCESS: ' + string, c='g')
