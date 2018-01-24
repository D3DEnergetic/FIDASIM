#!/usr/bin/env python
# -*- coding: utf-8 -*-

#+#FIDASIM Utilities
#+This file contains useful FIDASIM utilities
#+***
from __future__ import print_function
import os
from os.path import dirname
import subprocess
import platform
import numpy as np
import h5py

def get_fidasim_dir():
    """
    #+#get_fidasim_dir
    #+ Gets FIDASIM install directory
    #+***
    #+##Output Arguments
    #+     **directory**: FIDASIM install directory.
    #+##Example Usage
    #+```python
    #+>>> fida_dir = get_fidasim_dir()
    #+```
    """

    directory = dirname(dirname(dirname(dirname(os.path.abspath(__file__)))))

    return directory

def get_version(fidasim_dir):
    """
    #+#get_version
    #+ Gets FIDASIM version number from git.
    #+ Falls back to reading VERSION file when git is not available
    #+***
    #+##Input Arguments
    #+    **fidasim_dir**: FIDASIM install directory
    #+
    #+##Output Arguments
    #+     **version**: FIDAIM version number.
    #+
    #+##Example Usage
    #+```python
    #+>>> version = get_version(get_fidasim_dir())
    #+```
    """
    version = ''
    alt = False

    if platform.system() == 'Windows':
        alt = True
    else:
        # Location of .git folder
        git_dir = r'{}{}.git'.format(fidasim_dir, os.path.sep)

        # git is installed if git_file is a file
        proc = subprocess.Popen('command -v git', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        git_file = proc.communicate()[0].decode('utf-8')
        git_file = git_file.replace('\n', '')

        # Check that .git folder is present and git is installed
        if os.path.isfile(git_file) and os.path.isdir(git_dir):
            try:
                version = subprocess.check_output(['git', '--git-dir={}'.format(git_dir), 'describe', '--tags', '--always', '--dirty'])
                version = version.replace('\n', '')
            except:
                alt = True
        else:
            alt = True

    # If above didn't work, read version file
    if alt:
        # Git 'version' filepath
        ver_file = '{}{}VERSION'.format(fidasim_dir, os.path.sep)

        if os.path.isfile(ver_file):
            with open(ver_file) as f:
                version = f.read()

    return version

def rz_grid(rmin, rmax, nr, zmin, zmax, nz):
    """
    #+#rz_grid
    #+Creates interpolation grid
    #+***
    #+##Arguments
    #+    **rmin**: Minimum radius [cm]
    #+
    #+    **rmax**: Maximum radius [cm]
    #+
    #+    **nr**: Number of radii
    #+
    #+    **zmin**: Minimum Z value [cm]
    #+
    #+    **zmax**: Maximum Z value [cm]
    #+
    #+    **nz**: Number of z values
    #+
    #+##Return Value
    #+Interpolation grid dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> grid = rz_grid(0,200.0,200,-100,100,200)
    #+```
    """
    dr = (rmax - rmin) / (nr - 1)
    dz = (zmax - zmin) / (nz - 1)
    r = rmin + dr * np.arange(nr, dtype=np.float64)
    z = zmin + dz * np.arange(nz, dtype=np.float64)

    r2d = np.tile(r, (nz, 1)).T
    z2d = np.tile(z, (nr, 1))

    grid = {'r2d': r2d,
            'z2d': z2d,
            'r': r,
            'z': z,
            'nr': nr,
            'nz': nz}

    return grid

def colored(text, color): #, on_color=None, attrs=None):
    """
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
    """
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

def info(string):
    """
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
    """
    print(colored('INFO: ' + string, 'cyan'))

def warn(string):
    """
    #+#warn
    #+Print a warning message
    #+***
    #+##Arguments
    #+     **string**: message
    #+
    #+##Example Usage
    #+```python
    #+>>> warn("This may be a problem")
    #+```
    """
    print(colored('WARNING: ' + string, 'magenta'))

def error(string, halt=False):
    """
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
    """
    print(colored('ERROR: {}'.format(string), 'red'))

    if halt:
        raise Exception()

def success(string):
    """
    #+#success
    #+Print a success message
    #+***
    #+##Arguments
    #+     **string**: message
    #+
    #+##Example Usage
    #+```python
    #+>>> success("Yay!!!")
    #+```
    """
    print(colored('SUCCESS: ' + string, 'green'))

def beam_grid(nbi, rstart,
              nx=None, ny=None, nz=None, dv=8.0,
              length=100.0, width=80.0, height=80.0):
    """
    #+#beam_grid
    #+ Calculates settings for a grid that aligns with the neutral beam.
    #+***
    #+##Arguments
    #+    **nbi**: [Neutral beam geometry structure](|url|/page/03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure)
    #+
    #+    **rstart**: Radial start position of beam grid [cm]
    #+
    #+##Keyword Arguments
    #+    **dV**: Cell volume [\(cm^3\)]: Defaults to 8.0
    #+
    #+    **nx**: Number of cells in length: Default determined by `dV`
    #+
    #+    **ny**: Number of cells in width: Default determined by `dV`
    #+
    #+    **nz**: Number of cells in height: Default determined by `dV`
    #+
    #+    **length**: Length of grid along beam sightline. [cm]: Defaults to 100 cm
    #+
    #+    **width**: Width of grid [cm]: Defaults to 100 cm
    #+
    #+    **height**: Height of grid [cm]: Defaults 80 cm
    #+
    #+##Return Value
    #+    Structure containing beam grid settings suitable for the Namelist File
    #+
    #+##Example Usage
    #+```python
    #+>>> grid = beam_grid(nbi,200.0,nx=100,ny=50,nz=50,length=100,width=50,height=50)
    #+```
    """

    if width < nbi['widy']:
        warn("Grid width is smaller then the source width")
        print("width: {}".format(width))
        print("source width: {}".format(nbi['widy']))

    if height < nbi['widz']:
        warn("Grid height is smaller then the source height")
        print("height: {}".format(height))
        print("source height: {}".format(nbi['widz']))

    dv3 = dv ** (1. / 3.)

    if nx is None:
        nx = round(length / dv3)

    if ny is None:
        ny = round(width / dv3)

    if nz is None:
        nz = round(height / dv3)

    xmin = 0.
    xmax = length
    ymin = -width / 2.
    ymax = width / 2.
    zmin = -height / 2.
    zmax = height / 2.

    src = nbi['src']
    axis = nbi['axis'] / np.sqrt(np.sum(nbi['axis'] ** 2))
    pos = src + 100. * axis

    if np.sqrt(src[0] ** 2 + src[1] ** 2) < rstart:
        error("Source radius cannot be less then rstart", halt=True)

    dis = np.sqrt(np.sum((src - pos) ** 2.0))
    beta = np.arcsin((src[2] - pos[2]) / dis)
    alpha = np.arctan2((pos[1] - src[1]), (pos[0] - src[0]))
    gamma = 0.
    a = axis[0] ** 2 + axis[1] ** 2
    b = 2. * (src[0] * axis[0] + src[1] * axis[1])
    c = src[0] ** 2 + src[1] ** 2 - rstart ** 2
    t = (-b - np.sqrt(b ** 2 - 4. * a * c)) / (2. * a)
    origin = src + t * axis

    beam_grid = {'nx': nx,
                 'ny':  ny,
                 'nz':  nz,
                 'xmin': xmin,
                 'xmax': xmax,
                 'ymin': ymin,
                 'ymax': ymax,
                 'zmin': zmin,
                 'zmax': zmax,
                 'alpha': alpha,
                 'beta': beta,
                 'gamma': gamma,
                 'origin': origin}

    return beam_grid

def write_data(h5_obj, dic, desc, units, name=''):
    """
    #+#write_data
    #+ Write h5 datasets with attributes 'description' and 'units'
    #+***
    #+##Arguments
    #+     **h5_obj**: An h5 file or group object from h5py
    #+
    #+     **dic**: Dict of data to save as h5 datasets
    #+
    #+     **desc**: Dict with same keys as dic describing each item in dic
    #+
    #+     **units**: Dict with same keys as dic providing units of data in dic, doesn't have to be all keys of dic.
    #+
    #+##Keyword Arguments
    #+     **name**: Name/description of dic for clarity in raising errors
    #+
    #+##Example Usage
    #+```python
    #+>>> write_data(h5_obj, dic, desc, units)
    #+```
    """
    for key in dic:
        # Transpose data to match expected by Fortran and historically provided by IDL
        if isinstance(dic[key], np.ndarray):
            if dic[key].ndim >= 2:
                dic[key] = dic[key].T

        # Make strings of fixed length as required by Fortran.
        # See http://docs.h5py.org/en/latest/strings.html#fixed-length-ascii
        if isinstance(dic[key], str):
            dic[key] = np.string_(dic[key])

        # Create dataset
        ds = h5_obj.create_dataset(key, data = dic[key])

        # Add descrption attribute
        ds.attrs['description'] = desc[key]

        # Add units attribute (if present)
        if key in units:
            ds.attrs['units'] = units[key]
