#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from lib.python_prefida.info import info
from lib.python_prefida.error import error
from lib.python_prefida.check_dict_schema import check_dict_schema
from lib.python_prefida.success import success


def check_grid(grid):
    """
    ;+#check_grid
    ;+Checks if interpolation grid structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **grid**: Interpolation grid structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_grid, grid
    ;+```
    """
    err_status = 0
    info('Checking interpolation grid...')

    if 'nr' not in grid:
        error('"nr" is missing from the interpolation grid')
        err_status = 1
        error('Invalid interpolation grid. Exiting...', halt=True)

    if 'nz' not in grid:
        error('"nz" is missing from the interpolation grid')
        err_status = 1
        error('Invalid interpolation grid. Exiting...', halt=True)

    nr = grid['nr']
    nz = grid['nz']

    zero_int = {'dims': 0,
                'type': int}

    nrnz_doub = {'dims': [nr, nz],
                 'type': float}

    schema = {'nr': zero_int,
              'nz': zero_int,
              'r2d': nrnz_doub,
              'z2d': nrnz_doub,
              'r': {'dims': [nr],
                    'type': float},
              'z': {'dims': [nz],
                    'type': float}}

    err_status = check_dict_schema(schema, grid, desc="interpolation grid")
    if err_status == 1:
        error('Invalid interpolation grid. Exiting...', halt=True)

    if not np.array_equal(grid['r'], np.sort(grid['r'])):
        error('r is not in ascending order')
        err_status = 1

    if not np.array_equal(grid['z'], np.sort(grid['z'])):
        error('z is not in ascending order')
        err_status = 1

    if not np.array_equal(grid['r'], grid['r2d'][:, 0]):
        error('r2d is defined incorrectly. Expected r == r2d[:, 0]')
        err_status = 1

    if not np.array_equal(grid['z'], grid['z2d'][0, :]):
        error('z2d is defined incorrectly. Expected z == z2d[0, :]')
        err_status = 1

    if err_status != 0:
        error('Invalid interpolation grid. Exiting...', halt=True)
    else:
        success('Interpolation grid is valid')
