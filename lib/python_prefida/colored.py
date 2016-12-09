#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from termcolor import colored as col

#def colored(string, color):
#    return col(string, color=color)

import os

# The following is taken directly from termcolor package


def colored(text, color): #, on_color=None, attrs=None):
    """Colorize text.

    Available text colors:
        red, green, yellow, blue, magenta, cyan, white.

    Available text highlights:
        on_red, on_green, on_yellow, on_blue, on_magenta, on_cyan, on_white.

    Available attributes:
        bold, dark, underline, blink, reverse, concealed.

    Example:
        colored('Hello, World!', 'red', 'on_grey', ['blue', 'blink'])
        colored('Hello, World!', 'green')
    """
    COLORS = dict(list(zip(['grey',
                            'red',
                            'green',
                            'yellow',
                            'blue',
                            'magenta',
                            'cyan',
                            'white',],
                            list(range(30, 38)))))

    RESET = '\033[0m'

    if os.getenv('ANSI_COLORS_DISABLED') is None:
        fmt_str = '\033[%dm%s'
#        if color is not None:
        text = fmt_str % (COLORS[color], text)

#            if on_color is not None:
#                text = fmt_str % (HIGHLIGHTS[on_color], text)
#
#            if attrs is not None:
#                for attr in attrs:
#                    text = fmt_str % (ATTRIBUTES[attr], text)

        text += RESET

    return text
###############################################################################
if __name__ == "__main__":
    a = colored('test', 'green')

    print(a)
