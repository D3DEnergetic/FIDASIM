#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from python_prefida import info
from python_prefida import error
from python_prefida import check_dict_schema
from python_prefida import success


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

#    w = where("nr" eq strlowcase(TAG_NAMES(grid)),nw)
#    if nw eq 0:
    if 'nr' not in grid:
        error('"nr" is missing from the interpolation grid')
        err_status = 1
#        goto, GET_OUT
        error('Invalid interpolation grid. Exiting...', halt=True)


#    w = where("nz" eq strlowcase(TAG_NAMES(grid)),nw)
#    if nw eq 0:
    if 'nz' not in grid:
        error('"nz" is missing from the interpolation grid')
        err_status = 1
#        goto, GET_OUT
        error('Invalid interpolation grid. Exiting...', halt=True)

    nr = grid['nr']
    nz = grid['nz']
    zero_int = {'dims': 0, 'type': int}
    nrnz_doub = {'dims': [nr, nz], 'type': np.float64}
    schema = {'nr': zero_int, 'nz': zero_int,
              'r2d': nrnz_doub,
              'z2d': nrnz_doub,
              'r': {'dims': [nr], 'type': np.float64},
              'z': {'dims': [nz], 'type': np.float64}}

    err_status = check_dict_schema(schema, grid, desc="interpolation grid")
    if err_status == 1:
#        goto, GET_OUT
        error('Invalid interpolation grid. Exiting...', halt=True)

#    w = where((indgen(nr) eq sort(grid.r)) ne 1, nw)
#    if nw ne 0:
#        error,'r is not in ascending order'
#        err_status = 1
    if not np.array_equal(grid['r'], np.sort(grid['r'])):
        error('r is not in ascending order')
        err_status = 1

#    w = where((indgen(nz) eq sort(grid.z)) ne 1, nw)
#    if nw ne 0:
#        error,'z is not in ascending order'
#        err_status = 1
    if not np.array_equal(grid['z'], np.sort(grid['z'])):
        error('z is not in ascending order')
        err_status = 1

#    w = where((grid.r eq grid.r2d[*,0]) ne 1, nw)
#    if nw ne 0:
#        error,'r2d is defined incorrectly. Expected r == r2d[*,0]'
#        err_status = 1
    if not np.array_equal(grid['r'], grid['r2d'][:, 0]):
        error('r2d is defined incorrectly. Expected r == r2d[:, 0]')
        err_status = 1

#    w = where((grid.z eq grid.z2d[0,*]) ne 1, nw)
#    if nw ne 0:
#        error,'z2d is defined incorrectly. Expected z == z2d[0,*]'
#        err_status = 1
    if not np.array_equal(grid['z'], grid['z2d'][:, 0]):
        error('z2d is defined incorrectly. Expected z == z2d[:, 0]')
        err_status = 1

#    GET_OUT:
#    if err_status ne 0:
#        error,'Invalid interpolation grid. Exiting...',/halt
#     else:
#        success,'Interpolation grid is valid'

    if err_status != 0:
        error('Invalid interpolation grid. Exiting...', halt=True)
    else:
        success('Interpolation grid is valid')
