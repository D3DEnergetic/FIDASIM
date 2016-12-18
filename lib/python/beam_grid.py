#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
from lib.python_prefida.warn import warn
from lib.python_prefida.error import error


def beam_grid(nbi, rstart,
              nx=None, ny=None, nz=None, dv=None,
              length=None, width=None, height=None):
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
    #+    **width**: Width of grid [cm]: Defaults to 50 cm
    #+
    #+    **height**: Height of grid [cm]: Defaults to 50 cm
    #+
    #+##Return Value
    #+    Structure containing beam grid settings suitable for the Namelist File
    #+
    #+##Example Usage
    #+```python
    #+>>> grid = beam_grid(nbi,200.0,nx=100,ny=50,nz=50,length=100,width=50,height=50)
    #+```
    if length is None:
        length = 100.0      # cm

    if width is None:
        width = 80.0    # cm

    if width < nbi['widy']:
        warn("Grid width is smaller then the source width")
        print("width: {}".format(width))
        print("source width: {}".format(nbi['widy']))

    if height is None:
        height = 80.0   # cm

    if height < nbi['widz']:
        warn("Grid height is smaller then the source height")
        print("height: {}".format(height))
        print("source height: {}".format(nbi['widz']))

    if dv is None:
        dv = 8.0    # cm^3

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
