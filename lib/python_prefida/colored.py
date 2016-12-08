#!/usr/bin/env python
# -*- coding: utf-8 -*-

from termcolor import colored as col

def colored(string, color):

    return col(string, color=color)
###############################################################################
if __name__ == "__main__":
    a = colored('test', 'yellow')

    print(a)
