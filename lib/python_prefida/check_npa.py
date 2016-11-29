#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from lib import info
from lib import error
from lib import check_dict_schema
from lib import uvw_to_xyz
from lib import aabb_intersect
from lib import warn
from lib import success


def check_npa(inp, npa):
    """
    ;+#check_npa
    ;+Checks if NPA geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **npa**: NPA geometry structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_npa, inputs, npa
    ;+```
    """
    err_status = 0
    info('Checking NPA geometry...')

    npa_keys = npa.keys()

#    w = where("nchan" eq strlowcase(TAG_NAMES(npa)),nw)
#    if nw eq 0:
    if 'nchan' not in npa_keys:
        error('"nchan" is missing from the NPA geometry')
        err_status = 1
        error('Invalid NPA geometry. Exiting...', halt=True)

    nchan = npa['nchan']
    zero_string = {'dims': 0, 'type': str}
    zero_long = {'dims': 0, 'type': long}
    schema = {'data_source': zero_string,
              'nchan': zero_long,
              'system': zero_string,
              'id': {'dims': [nchan], 'type': str},
              'a_shape': {'dims': [nchan], 'type': int},
              'd_shape': {'dims': [nchan], 'type': int},
              'a_tedge': {'dims': [3, nchan], 'type': np.float64},
              'a_redge': {'dims': [3, nchan], 'type': np.float64},
              'a_cent': {'dims': [3, nchan], 'type': np.float64},
              'd_tedge': {'dims': [3, nchan], 'type': np.float64},
              'd_redge': {'dims': [3, nchan], 'type': np.float64},
              'd_cent': {'dims': [3, nchan], 'type': np.float64},
              'radius': {'dims': [nchan], 'type': np.float64}}

    err_status = check_dict_schema(schema, npa, desc="NPA geometry")
    if err_status == 1:
        error('Invalid NPA geometry. Exiting...', halt=True)

    # Check detector/aperture shape
    w = (npa['d_shape'] > 2) or (npa['d_shape'] == 0)
    nw = len(npa['d_shape'][w])
    if nw != 0:
#    if npa['d_shape'] not in [1, 2]:
        error('Invalid detector shape. Expected 1 (rectagular) or 2 (circular)')
        print 'Invalid indices: {}'.format(np.arange(len(npa['d_shape']))[w])
        err_status = 1

    w = (npa['a_shape'] > 2) or (npa['a_shape'] == 0)
    nw = len(npa['a_shape'])
    if nw != 0:
        error('Invalid aperture shape. Expected 1 (rectagular) or 2 (circular)')
        print 'Invalid indices: {}'.format(np.arange(len(npa['a_shape']))[w])
        err_status = 1

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
#    print f='(i3," out of ",i3," channels crossed the beam grid")',nw,nchan
    print '{} out of {} channels crossed the beam grid'.format(nw, nchan)
    if nw == 0:
        error('No channels intersect the beam grid')
        err_status = 1

    if nww > 0:
        warn('Some channels did not intersect the beam grid')
        print 'Number missed: {}'.format(nww)
        print 'Missed channels:'
        print '    {}'.format(npa['id'][ww])

#    GET_OUT:
    if err_status != 0:
        error('Invalid NPA geometry. Exiting...', halt=True)
    else:
        success('NPA geometry is valid')
