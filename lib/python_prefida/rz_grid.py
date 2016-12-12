#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def rz_grid(rmin, rmax, nr, zmin, zmax, nz):
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
