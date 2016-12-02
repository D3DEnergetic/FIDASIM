#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from python_prefida import info
from python_prefida import check_dict_schema
from python_prefida import error
from python_prefida import uvw_to_xyz
from python_prefida import aabb_intersect
from python_prefida import success
import numpy as np


def check_beam(inp, nbi):
    """
    #check_beam
    Checks if neutral beam geometry structure is valid
    ***
    ##Input Arguments
         **inputs**: input structure

         **nbi**: neutral beam geometry structure

    ##Example Usage
    ```idl
    IDL> check_beam, inputs, nbi
    ```
    """
    err_status = 0
    info('Checking beam geometry...')

    na = nbi['naperture']
    zero_string = {'dims': 0, 'type': str}
    zero_int = {'dims': 0, 'type': int}
    zero_double = {'dims': 0, 'type': np.float64}
    three_double = {'dims': [3], 'type': np.float64}
    na_double = {'dims': [na], 'type': np.float64}
    na_int = {'dims': [na], 'type': int}
    schema = {'data_source': zero_string,
              'name': zero_string, 'shape': zero_int,
              'src': three_double, ' axis': three_double,
              'divy': three_double, 'divz': three_double,
              'focy': zero_double, ' focz': zero_double,
              'widz': zero_double, ' widy': zero_double,
              'naperture': zero_int, 'ashape': na_int,
              'awidy': na_double, 'awidz': na_double,
              'aoffy': na_double, 'aoffz': na_double,
              'adist': na_double}

    err_status = check_dict_schema(schema, nbi, desc="beam geometry")
    if err_status:
#        goto, GET_OUT
        error('Invalid beam geometry. Exiting...', halt=True)

    if np.abs(np.sum(nbi['axis'] ** 2.) - 1.) > 1e-5:
        error('Invalid source axis. Expected norm(axis) == 1')
        err_status = 1

    if nbi['focz'] <= 0.0:
        error('focz cannot be in the range (-Inf,0.0]')
        err_status = 1

    if nbi['focy'] <= 0.0:
        error('focy cannot be in the range (-Inf,0.0]')
        err_status = 1

    if nbi['shape'] not in [1, 2]:
        error('Invalid source shape. Expected 1 (rectagular) or 2 (circular)')
        err_status = 1

    if nbi['widz'] < 0.:
        error('Invalid widz. Expected widz > 0')
        err_status = 1

    if nbi['widy'] < 0:
        error('Invalid widy. Expected widy > 0')
        err_status = 1

#    w = where(nbi.ashape > 2 or nbi.ashape eq 0,nw)
    if nbi['ashape'] not in [1, 2]:
        error('Invalid aperture shape. Expected 1 (rectangular) or 2 (circular)')
        err_status = 1

    w = nbi['awidy'] < 0
    nw = len(nbi['awidy'][w])
    if nw > 0:
        error('Invalid awidy. Expected awidy >= 0.0')
        err_status = 1

#    w = where(nbi.awidz < 0, nw)
    w = nbi['awidz'] < 0
    nw = len(nbi['awidz'][w])
    if nw > 0:
        error('Invalid awidz. Expected awidz >= 0.0')
        err_status = 1

    origin = inp['origin']
    uvw_src = nbi['src']
    uvw_axis = nbi['axis']
    uvw_pos = uvw_src + nbi['adist'][0] * uvw_axis

    xyz_src = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_src, origin=origin)
    xyz_axis = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_axis)
    xyz_pos = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_pos, origin=origin)
    xyz_center = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], [0., 0., 0.], origin=origin)

    dis = np.sqrt(np.sum((xyz_src - xyz_pos) ** 2.))
    BETA = np.float64(np.arcsin((xyz_src[2] - xyz_pos[2]) / dis))
    ALPHA = np.float64(np.arctan((xyz_pos[1] - xyz_src[1]), (xyz_pos[0] - xyz_src[0])))

#    print('Beam injection start point in machine coordinates'
#    print( f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_src
#    print('First aperture position in machine coordinates'
#    print( f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_pos
#    print('Machine center in beam grid coordinates'
#    print( f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_center
#    print('Beam injection start point in beam grid coordinates'
#    print( f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_src
#    print('First aperture position in beam grid coordinates'
#    print( f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_pos

    print('Beam injection start point in machine coordinates')
    print(uvw_src)
    print('First aperture position in machine coordinates')
    print(uvw_pos)
    print('Machine center in beam grid coordinates')
    print(xyz_center)
    print('Beam injection start point in beam grid coordinates')
    print(xyz_src)
    print('First aperture position in beam grid coordinates')
    print(xyz_pos)

    print('Beam grid rotation angles that would align it with the beam centerline')
    print('{} deg.'.format(ALPHA / np.pi * 180.))  # ,FORMAT='("    alpha = ",F14.10,"°")'
    print('{} deg.'.format(BETA / np.pi * 180.))  # ,FORMAT='("    beta = ",F14.10,"°")'

    # Calculate grid center rc and sides length dr
    dr = np.array([inp['xmax'] - inp['xmin'], inp['ymax'] - inp['ymin'], inp['zmax'] - inp['zmin']], dtype=np.float64)
    rc = np.array([inp['xmin'], inp['ymin'], inp['zmin']], dtype=np.float64) + 0.5 * dr

    # Check if beam centerline intersects beam grid
#    aabb_intersect,rc,dr,xyz_src,xyz_axis,length,r_enter,r_exit
    length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_src, xyz_axis)

    print('Beam centerline - grid intersection length')
#    print(f='("    length = ",F8.3)',length
    print(length)
    if length <= 10.0:
        error('Beam centerline does not intersect grid')
        err_status = 1

#    GET_OUT:
    if err_status != 0:
        error('Invalid beam geometry. Exiting...', halt=True)
    else:
        success('Beam geometry is valid')
