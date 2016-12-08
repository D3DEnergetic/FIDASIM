#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from lib.python_prefida.info import info
from lib.python_prefida.error import error
from lib.python_prefida.warn import warn
from lib.python_prefida.success import success
from lib.python_prefida.check_dict_schema import check_dict_schema


def check_fields(inp, grid, fields):
    #+#check_fields
    #+Checks if electromagnetic fields structure is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: Input structure
    #+
    #+     **grid**: Interpolation grid structure
    #+
    #+     **fields**: Electromagnetic fields structure
    #+
    #+##Example Usage
    #+```idl
    #+IDL> check_fields, inputs, grid, fields
    #+```
    err = False
    info('Checking electromagnetic fields...')

    nr = grid['nr']
    nz = grid['nz']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_double = {'dims': 0,
                   'type': [float, np.float64]}

    nrnz_double = {'dims': [nr, nz],
                   'type': [float, np.float64]}

    nrnz_int = {'dims': [nr, nz],
                'type': [int, np.int32]}

    schema = {'time': zero_double,
              'br': nrnz_double,
              'bt': nrnz_double,
              'bz': nrnz_double,
              'er': nrnz_double,
              'et': nrnz_double,
              'ez': nrnz_double,
              'mask': nrnz_int,
              'data_source': zero_string}

    err = check_dict_schema(schema, fields, desc="electromagnetic fields")
    if err:
        error('Invalid electromagnetic fields. Exiting...', halt=True)

    if fields['data_source'] == '':
        error('Invalid data source. An empty string is not a data source.')
        err = True

    if np.abs(fields['time'] - inp['time']) > 0.02:
        warn('Electromagnetic fields time and input time do not match')
        print('Input time: {}'.format(inp['time']))
        print('Electromagnetic fields time: {}'.format(fields['time']))

    # Add grid elements to fields dict
    for key in grid:
        fields[key] = grid[key]

    if err:
        error('Invalid electromagnetic fields. Exiting...', halt=True)
    else:
        success('Electromagnetic fields are valid')

    return fields
