#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
from lib.python_prefida import info
from lib.python_prefida import error
from lib.python_prefida import check_dict_schema
from lib.python_prefida import uvw_to_xyz
from lib.python_prefida import aabb_intersect
from lib.python_prefida import warn
from lib.python_prefida import success


def check_npa(inp, npa):
    #+#check_npa
    #+Checks if NPA geometry dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: input dictionary
    #+
    #+     **npa**: NPA geometry dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> check_npa(inputs, npa)
    #+```
    err = False
    info('Checking NPA geometry...')

    npa_keys = npa.keys()

    if 'nchan' not in npa_keys:
        error('"nchan" is missing from the NPA geometry')
        err = True
        error('Invalid NPA geometry. Exiting...', halt=True)

    nchan = npa['nchan']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_long = {'dims': 0,
                 'type': [int, np.int32, np.int64]}

    three_float = {'dims': [3, nchan],
                   'type': [float, np.float64]}

    nchan_int = {'dims': [nchan],
                 'type': [int, np.int32, np.int64]}

    schema = {'data_source': zero_string,
              'nchan': zero_long,
              'system': zero_string,
              'id': {'dims': [nchan],
                     'type': [str]},
              'a_shape': nchan_int,
              'd_shape': nchan_int,
              'a_tedge': three_float,
              'a_redge': three_float,
              'a_cent': three_float,
              'd_tedge': three_float,
              'd_redge': three_float,
              'd_cent': three_float,
              'radius': three_float}

    err = check_dict_schema(schema, npa, desc="NPA geometry")
    if err:
        error('Invalid NPA geometry. Exiting...', halt=True)

    # Check detector/aperture shape
    w = (npa['d_shape'] > 2) or (npa['d_shape'] == 0)
    nw = len(npa['d_shape'][w])
    if nw != 0:
        error('Invalid detector shape. Expected 1 (rectagular) or 2 (circular)')
        print('Invalid indices: {}'.format(np.arange(len(npa['d_shape']))[w]))
        err = True

    w = (npa['a_shape'] > 2) or (npa['a_shape'] == 0)
    nw = len(npa['a_shape'])
    if nw != 0:
        error('Invalid aperture shape. Expected 1 (rectagular) or 2 (circular)')
        print('Invalid indices: {}'.format(np.arange(len(npa['a_shape']))[w]))
        err = True

    # Calculate grid center rc and sides length dr
    dr = [inp['xmax'] - inp['xmin'], inp['ymax'] - inp['ymin'], inp['zmax'] - inp['zmin']]
    rc = [inp['xmin'], inp['ymin'], inp['zmin']] + 0.5 * dr
    err_arr = np.zeros(nchan, dtype=int)
    for i in range(nchan):
        uvw_det = npa['d_cent'][:, i]
        d_e1 = npa['d_redge'][:, i] - uvw_det
        d_e2 = npa['d_tedge'][:, i] - uvw_det

        uvw_aper = npa['a_cent'][:, i]
        a_e1 = npa['a_redge'][:, i] - uvw_aper
        a_e2 = npa['a_tedge'][:, i] - uvw_aper

        #Rotate chords into beam grid coordinates
        xyz_aper = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_aper, origin=inp['origin'])
        xyz_det = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_det, origin=inp['origin'])
        xyz_dir = xyz_aper - xyz_det
        xyz_dir = xyz_dir / np.sqrt(np.sum(xyz_dir * xyz_dir))

        # Check if npa chord intersects beam grid
        length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_det, xyz_dir)
        if length <= 0.:
            err_arr[i] = 1

        # Check that the detector and aperture point in the same direction
        d_e3 = np.cross(d_e1, d_e2)
        a_e3 = np.cross(a_e1, a_e2)
        dp = np.sum(d_e3 * a_e3)
        if dp <= 0.:
            warn('The dot product of the detector and aperture plane normal vectors is negative. The NPA definition may be incorrect.')

    w = (err_arr == 0)
    nw = err_arr[w].size
    ww = (err_arr != 0)
    nww = err_arr[ww].size
    print('{} out of {} channels crossed the beam grid'.format(nw, nchan))
    if nw == 0:
        error('No channels intersect the beam grid')
        err = True

    if nww > 0:
        warn('Some channels did not intersect the beam grid')
        print('Number missed: {}'.format(nww))
        print('Missed channels:')
        print('    {}'.format(npa['id'][ww]))

    if err:
        error('Invalid NPA geometry. Exiting...', halt=True)
    else:
        success('NPA geometry is valid')
