#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import IPython.utils.path
import os
from lib.python_prefida.info import info
from lib.python_prefida.check_dict_schema import check_dict_schema
from lib.python_prefida.error import error
from lib.python_prefida.success import success


def check_inputs(inputs):
    #+#check_inputs
    #+Checks if input structure is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: input structure
    #+
    #+##Example Usage
    #+```idl
    #+IDL> check_inputs, inputs
    #+```
    info('Checking simulation settings...')
    err = False

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_int = {'dims': 0,
                'type': [int]}

    zero_long = {'dims': 0,
                 'type': [int]}

    zero_double = {'dims': 0,
                   'type': [float, np.float64]}

    three_double = {'dims': [3],
                    'type': [float, np.float64]}

    schema = {'comment': zero_string,
              'shot': zero_long,
              'time': zero_double,
              'runid': zero_string,
              'device': zero_string,
              'tables_file': zero_string,
              'result_dir': zero_string,
              'nlambda': zero_int,
              'lambdamin': zero_double,
              'lambdamax': zero_double,
              'nx': zero_int,
              'ny': zero_int,
              'nz': zero_int,
              'alpha': zero_double,
              'beta': zero_double,
              'gamma': zero_double,
              'origin': three_double,
              'xmin': zero_double,
              'xmax': zero_double,
              'ymin': zero_double,
              'ymax': zero_double,
              'zmin': zero_double,
              'zmax': zero_double,
              'ab': zero_double,
              'ai': zero_double,
              'current_fractions': three_double,
              'pinj': zero_double,
              'einj': zero_double,
              'impurity_charge': zero_int,
              'n_fida': zero_long,
              'n_nbi': zero_long,
              'n_dcx': zero_long,
              'n_npa': zero_long,
              'n_halo': zero_long,
              'n_birth': zero_long,
              'ne_wght': zero_int,
              'np_wght': zero_int,
              'nphi_wght': zero_int,
              'emax_wght': zero_double,
              'nlambda_wght': zero_int,
              'lambdamin_wght': zero_double,
              'lambdamax_wght': zero_double,
              'calc_npa': zero_int,
              'calc_fida': zero_int,
              'calc_bes': zero_int,
              'calc_brems': zero_int,
              'calc_birth': zero_int,
              'calc_fida_wght': zero_int,
              'calc_npa_wght': zero_int,
              'dump_dcx': zero_int}

    err = check_dict_schema(schema, inputs, desc="simulation settings")
    if err:
        error('Invalid simulation settings. Exiting...', halt=True)

    # Normalize File Paths
    inputs['result_dir'] = IPython.utils.path.expand_path(inputs['result_dir'])

    if (inputs['alpha'] > 2. * np.pi) or (inputs['beta'] > 2. * np.pi) or (inputs['gamma'] > 2. * np.pi):
        error('Angles must be in radians')
        err = True

    if inputs['lambdamin'] >= inputs['lambdamax']:
        error('Invalid wavelength range. Expected lambdamin < lamdbdamax')
        err = True

    if inputs['lambdamin_wght'] >= inputs['lambdamax_wght']:
        error('Invalid wavelength range. Expected lambdamin_wght < lamdbdamax_wght')
        err = True

    if inputs['xmin'] >= inputs['xmax']:
        error('Invalid x range. Expected xmin < xmax')
        err = True

    if inputs['ymin'] >= inputs['ymax']:
        error('Invalid y range. Expected ymin < ymax')
        err = True

    if inputs['zmin'] >= inputs['zmax']:
        error('Invalid z range. Expected zmin < zmax')
        err = True

    if (inputs['pinj'] <= 0.) or (inputs['einj'] <= 0.0):
        error('The selected source is not on')
        print('einj = {}'.format(inputs['einj']))
        print('pinj = {}'.format(inputs['pinj']))
        err = True

    if np.abs(np.sum(inputs['current_fractions']) - 1.0) > 1e-3:
        error('current_fractions do not sum to 1.0')
        print('sum(current_fractions) = {}'.format(np.sum(inputs['current_fractions'])))
        err = True

    if inputs['impurity_charge'] <= 1:
        error('Invalid impurity charge. Expected impurity charge > 1')
        err = True

    ps = os.path.sep
    input_file = inputs['result_dir'] + ps + inputs['runid'] + '_inputs.dat'
    equilibrium_file = inputs['result_dir'] + ps + inputs['runid'] + '_equilibrium.h5'
    geometry_file = inputs['result_dir'] + ps + inputs['runid'] + '_geometry.h5'
    distribution_file = inputs['result_dir'] + ps + inputs['runid'] + '_distribution.h5'
    neutrals_file = inputs['result_dir'] + ps + inputs['runid'] + '_neutrals.h5'

    inputs['input_file'] = input_file
    inputs['equilibrium_file'] = equilibrium_file
    inputs['geometry_file'] = geometry_file
    inputs['distribution_file'] = distribution_file
    inputs['load_neutrals'] = 0
    inputs['no_flr'] = 0
    inputs['verbose'] = 1
    inputs['neutrals_file'] = neutrals_file

    if err:
        error('Invalid simulation settings. Exiting...', halt=True)
    else:
        success('Simulation settings are valid')

    return inputs
