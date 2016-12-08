#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib.python_prefida.info import info
from lib.python_prefida.check_dict_schema import check_dict_schema
from lib.python_prefida.uvw_to_xyz import uvw_to_xyz
from lib.python_prefida.error import error
from lib.python_prefida.warn import warn
from lib.python_prefida.success import success
from lib.python_prefida.aabb_intersect import aabb_intersect
import numpy as np


def check_spec(inputs, chords):
    """
    #check_spec
    Check if spectral geometry structure is valid
    ***
    ##Input Arguments
         **inputs**: input structure

         **chords**: spectral geometry structure

    ##Example Usage
    ```idl
    IDL> check_spec, inputs, chords
    ```
    """
    err_status = False
    info('Checking FIDA/BES inputs...')

    chords_keys = list(chords.keys())

#    w = where("nchan" == strlowcase(TAG_NAMES(chords)),nw)
#    if nw == 0:
    if 'nchan' not in chords_keys:
        error('"nchan" is missing from the FIDA/BES geometry')
        err_status = True
        error('Invalid FIDA/BES geometry. Exiting...', halt=True)

    nchan = chords['nchan']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_long = {'dims': 0,
                 'type': [int, np.int32]}

    nchan_double = {'dims': [nchan],
                    'type': [float, np.float64]}

    nchan_string = {'dims': [nchan],
                    'type': [str, np.str_, np.bytes_]}

    schema = {'data_source': zero_string,
              'nchan': zero_long,
              'system': zero_string,
              'id': nchan_string,
              'lens': {'dims': [3, nchan],
                       'type': [float, np.float64]},
              'axis': {'dims': [3, nchan],
                       'type': [float, np.float64]},
              'sigma_pi': nchan_double,
              'spot_size': nchan_double,
              'radius': nchan_double}

    err_status = check_dict_schema(schema, chords, desc="FIDA/BES geometry")
    if err_status:
        error('Invalid FIDA/BES geometry. Exiting...', halt=True)

#    err_arr = np.zeros(nchan, dtype=int)
    cross_arr = np.zeros(nchan, dtype=int)
    uvw_lens = chords['lens']
    uvw_axis = chords['axis']

    #ROTATE CHORDS INTO BEAM GRID COORDINATES
    xyz_lens = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_lens, origin=inputs['origin'])
    xyz_axis = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_axis)

    # Calculate grid center rc and sides length dr
    dr = np.array([inputs['xmax'] - inputs['xmin'], inputs['ymax'] - inputs['ymin'], inputs['zmax'] - inputs['zmin']], dtype=float)
    rc = np.array([inputs['xmin'], inputs['ymin'], inputs['zmin']], dtype=float) + 0.5 * dr

    for i in range(nchan):
#        chan_str = strcompress(string(i),/remove_all)
        if np.abs(np.sum(uvw_axis[:, i] ** 2.) - 1.) > 1e-5:
            error('Invalid optical axis for chord "' + chords['id'][i] + '". Expected norm(axis) == 1')
            print(np.sum(uvw_axis[:, i] ** 2.) - 1.)
#            err_arr[i] = 1

        # Check if viewing chord intersects beam grid
        length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_lens[:, i], xyz_axis[:, i])
        if length <= 0.0:
            warn('Chord "{}" does not cross the beam grid'.format(chords['id'][i]))
#            warn('Chord "' + chords['id'][i] + '" does not cross the beam grid')
            cross_arr[i] = 1

#    w = where(cross_arr == 0.0, nw, complement=ww, ncomplement=nww)
#    print f='(i3," out of ",i3," chords crossed the beam grid")',nw,nchan
    w = (cross_arr == 0)
    nw = cross_arr[w].size
    print('{} out of {} chords crossed the beam grid'.format(nw, nchan))
    if nw == 0:
#    if 1 not in cross_arr:
        error('No channels intersect the beam grid')
        err_status = True

#    GET_OUT:
    if err_status:
        error('Invalid FIDA/BES geometry. Exiting...', halt=True)
    else:
        success('FIDA/BES geometry is valid')
