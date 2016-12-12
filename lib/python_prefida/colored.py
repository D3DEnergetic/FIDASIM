#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os


def colored(text, color): #, on_color=None, attrs=None):
    #+#colored
    #+ Return text string formatting for color in terminal
    #+***
    #+##Input Arguments
    #+     **text**: String to be colored
    #+
    #+     **color**: Desired color of string. Red, green, yellow, blue, magenta, cyan, or white.
    #+
    #+##Output Arguments
    #+     **text**: Text formated to have "color" in terminal.
    #+##Example Usage
    #+```python
    #+>>> text = colored("Text to be red", 'red')
    #+>>> print(text)
    #+```
    # Copyright (c) 2008-2011 Volvox Development Team
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    # THE SOFTWARE.
    #
    # Author: Konstantin Lepa <konstantin.lepa@gmail.com>
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

        text = fmt_str % (COLORS[color], text)

        text += RESET

    return text
