#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from lib.python_prefida.success import success
from lib.python_prefida.info import info
from lib.python_prefida.check_dict_schema import check_dict_schema
from lib.python_prefida.error import error
from lib.python_prefida.warn import warn


def check_plasma(inp, grid, plasma):
    """
    #check_plasma
    Checks if plasma paramters structure is valid
    ***
    ##Input Arguments
         **inputs**: Input structure

         **grid**: Interpolation grid structure

         **plasma**: Plasma parameters structure

    ##Example Usage
    ```idl
    IDL> check_plasma, inputs, grid, plasma
    ```
    """
    err_status = 0
    info('Checking plasma parameters...')

    nr = grid['nr']
    nz = grid['nz']

    zero_string = {'dims': 0,
                   'type': str}

    zero_double = {'dims': 0,
                   'type': float}

    nrnz_double = {'dims': [nr, nz],
                   'type': float}

    nrnz_int = {'dims': [nr, nz],
                'type': int}

    schema = {'time': zero_double,
              'vr': nrnz_double,
              'vt': nrnz_double,
              'vz': nrnz_double,
              'dene': nrnz_double,
              'ti': nrnz_double,
              'te': nrnz_double,
              'zeff': nrnz_double,
              'mask': nrnz_int,
              'data_source': zero_string}

    err_status = check_dict_schema(schema, plasma, desc="plasma parameters")
    if err_status == 1:
        error('Invalid plasma parameters. Exiting...', halt=True)

    if plasma['data_source'] == '':
        error, 'Invalid data source. An empty string is not a data source.'
        err_status = 1

    # Electron density
    w = (plasma['dene'] < 0.)
    plasma['dene'][w] = 0.  # [1/cm^3]

    # Zeff
    w = (plasma['zeff'] < 1.)
    # plasma['zeff'] = plasma['zeff'] > 1.0
    plasma['zeff'][w] = 1.

    # Electron temperature
    #plasma['te'] = plasma['te'] > 0.  # [keV]
    w = (plasma['te'] < 0.)
    plasma['te'][w] = 0.

    # Ion temperature
    #plasma['ti'] = plasma['ti'] > 0.  # [keV]
    w = (plasma['ti'] < 0.)
    plasma['ti'][w] = 0.

    if (np.abs(plasma['time'] - inp['time']) > 0.02):
        warn('Plasma time and input time do not match')
        print('Input time: ', inp['time'])
        print('Plasma time: ', plasma['time'])

    #plasma = create_struct(plasma, grid)
    plasma['grid'] = grid

    if err_status != 0:
        error('Invalid plasma parameters. Exiting...', halt=True)
    else:
        success('Plasma parameters are valid')

    return plasma
