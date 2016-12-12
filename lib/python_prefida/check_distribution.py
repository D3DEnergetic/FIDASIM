#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib.python_prefida.info import info
from lib.python_prefida.check_dict_schema import check_dict_schema
from lib.python_prefida.warn import warn
from lib.python_prefida.error import error
from lib.python_prefida.success import success
import numpy as np


def check_distribution(inputs, grid, dist):
    #+#check_distribution
    #+Checks if distribution dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: Input dictionary
    #+
    #+     **grid**: Interpolation grid dictionary
    #+
    #+     **dist**: Fast-ion distribution dictionary
    #+
    #+##Output Arguments
    #+     **dist**: Updated dist dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> dist = check_distribution(inputs, grid, dist)
    #+```
    err = False
    info('Checking fast-ion distribution...')

    dist_keys = list(dist.keys())

    if 'type' not in dist_keys:
        error('"type" is missing from the fast-ion distribution')
        error('Invalid fast-ion distribution. Exiting...', halt=True)

    dist_type = dist['type']

    if dist_type == 1:
        print('Using a Guiding Center Fast-ion Density Function')
        if 'nenergy' not in dist_keys:
            error('"nenergy" is missing from the fast-ion distribution')
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        if 'npitch' not in dist_keys:
            error('"npitch" is missing from the fast-ion distribution')
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        npitch = dist['npitch']
        nen = dist['nenergy']
        nr = grid['nr']
        nz = grid['nz']

        zero_string = {'dims': 0,
                       'type': [str]}

        zero_int = {'dims': 0,
                    'type': [int, np.int32, np.int64]}

        zero_double = {'dims': 0,
                       'type': [float, np.float64]}

        nrnz_double = {'dims': [nr, nz],
                       'type': [float, np.float64]}

        schema = {'type': zero_int,
                  'nenergy': zero_int,
                  'npitch': zero_int,
                  'energy': {'dims': [nen],
                             'type': [float, np.float64]},
                  'pitch': {'dims': [npitch],
                            'type': [float, np.float64]},
                  'denf': nrnz_double,
                  'f': {'dims': [nen, npitch, nr, nz],
                        'type': [float, np.float64]},
                  'time': zero_double,
                  'data_source': zero_string}

        err = check_dict_schema(schema, dist, desc="fast-ion distribution")
        if err:
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        dist['grid'] = grid
    elif dist_type == 2:
        print('Using Guiding Center Monte Carlo fast-ion distribution')
        if 'nparticle' not in dist_keys:
            error('"nparticle" is missing from the fast-ion distribution')
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        npart = dist['nparticle']

        zero_int = {'dims': 0,
                    'type': [int, np.int32, np.int64]}

        zero_long = {'dims': 0,
                     'type': [int, np.int32, np.int64]}

        zero_string = {'dims': 0,
                       'type': [str]}

        zero_double = {'dims': 0,
                       'type':[float,  np.float64]}

        npart_double = {'dims': [npart],
                        'type': [float, np.float64]}

        npart_int = {'dims': [npart],
                     'type': [int, np.int32, np.int64]}

        schema = {'type': zero_int,
                  'nparticle': zero_long,
                  'nclass': zero_int,
                  'time': zero_double,
                  'energy': npart_double,
                  'pitch': npart_double,
                  'r': npart_double,
                  'z': npart_double,
                  'weight': npart_double,
                  'class': npart_int,
                  'data_source': zero_string}

        err = check_dict_schema(schema, dist, desc="fast-ion distribution")
        if err:
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        print('Number of MC particles: {}'.format(npart))
    elif dist_type == 3:
        print('Using Full Orbit Monte Carlo fast-ion distribution')
        if 'nparticle' not in dist_keys:
            error('"nparticle" is missing from the fast-ion distribution')
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        npart = dist['nparticle']

        zero_int = {'dims': 0,
                    'type': [int, np.int32, np.int64]}

        zero_long = {'dims': 0,
                     'type': [int, np.int32, np.int64]}

        zero_string = {'dims': 0,
                       'type': [str]}

        zero_double = {'dims': 0,
                       'type': [float, np.float64]}

        npart_double = {'dims': [npart],
                        'type': [float, np.float64]}

        npart_int = {'dims': [npart],
                     'type': [int, np.int32, np.int64]}

        schema = {'type': zero_int,
                  'nparticle': zero_long,
                  'nclass': zero_int,
                  'time': zero_double,
                  'vr': npart_double,
                  'vt': npart_double,
                  'vz': npart_double,
                  'r': npart_double,
                  'z': npart_double,
                  'weight': npart_double,
                  'class': npart_int,
                  'data_source': zero_string}

        err = check_dict_schema(schema, dist, desc="fast-ion distribution")
        if err:
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        print('Number of MC particles: {}'.format(npart))
    else:
        error('Invalid distribution type. Expected ' +
              '1 (Guiding Center Density Function), ' +
              '2 (Guiding Center Monte Carlo), or ' +
              '3 (Full Orbit Monte Carlo)')
        error('Invalid fast-ion distribution. Exiting...', halt=True)

    if dist['data_source'] == '':
        error('Invalid data source. An empty string is not a data source.')
        err = True

    if np.abs(dist['time'] - inputs['time']) > 0.02:
        warn('Distribution time and input time do not match')
        print('Input time: {}'.format(inputs['time']))
        print('Distribution time: {}'.format(dist['time']))

    if err:
        error('Invalid fast-ion distribution. Exiting...', halt=True)
    else:
        success('Fast-ion distribution is valid')

    return dist
