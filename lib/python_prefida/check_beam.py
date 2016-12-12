#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.info import info
from lib.python_prefida.check_dict_schema import check_dict_schema
from lib.python_prefida.error import error
from lib.python_prefida.uvw_to_xyz import uvw_to_xyz
from lib.python_prefida.aabb_intersect import aabb_intersect
from lib.python_prefida.success import success
import numpy as np


def check_beam(inputs, nbi):
    #+#check_beam
    #+Checks if neutral beam geometry dictionary is valid. Converts lists to numpy ndarrays
    #+***
    #+##Input Arguments
    #+     **inputs**: input dictionary
    #+
    #+     **nbi**: neutral beam geometry dictionary
    #+
    #+##Output Arguments
    #+     **nbi**: Updated nbi dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> nbi = check_beam(inputs, nbi)
    #+```
    err = False
    info('Checking beam geometry...')

    na = nbi['naperture']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_int = {'dims': 0,
                'type': [int, np.int32, np.int64]}

    zero_double = {'dims': 0,
                   'type': [float, np.float64]}

    three_double = {'dims': [3],
                    'type': [float, np.float64]}

    na_double = {'dims': [na],
                 'type': [float, np.float64]}

    na_int = {'dims': [na],
              'type': [int, np.int32, np.int64]}

    schema = {'data_source': zero_string,
              'name': zero_string,
              'shape': zero_int,
              'src': three_double,
              'axis': three_double,
              'divy': three_double,
              'divz': three_double,
              'focy': zero_double,
              'focz': zero_double,
              'widz': zero_double,
              'widy': zero_double}

    # Add to schema if aperatures are present
    if nbi['naperture'] > 0:
        schema['naperture'] = zero_int
        schema['ashape'] = na_int
        schema['awidy'] = na_double
        schema['awidz'] = na_double
        schema['aoffy'] = na_double
        schema['aoffz'] = na_double
        schema['adist'] = na_double

    # Convert to np arrays for indexing
    nbi['ashape'] = np.array(nbi['ashape'], dtype=int)
    nbi['awidy'] = np.array(nbi['awidy'], dtype=np.float64)
    nbi['awidz'] = np.array(nbi['awidz'], dtype=np.float64)
    nbi['aoffy'] = np.array(nbi['aoffy'], dtype=np.float64)
    nbi['aoffz'] = np.array(nbi['aoffz'], dtype=np.float64)
    nbi['adist'] = np.array(nbi['adist'], dtype=np.float64)

    err = check_dict_schema(schema, nbi, desc="beam geometry")
    if err:
        error('Invalid beam geometry. Exiting...', halt=True)

    if np.abs(np.sum(nbi['axis'] ** 2.) - 1.) > 1e-5:
        error('Invalid source axis. Expected norm(axis) == 1')
        err = True

    if nbi['focz'] <= 0.0:
        error('focz cannot be in the range (-Inf,0.0]')
        err = True

    if nbi['focy'] <= 0.0:
        error('focy cannot be in the range (-Inf,0.0]')
        err = True

    if nbi['shape'] not in [1, 2]:
        error('Invalid source shape. Expected 1 (rectagular) or 2 (circular)')
        err = True

    if nbi['widz'] < 0.:
        error('Invalid widz. Expected widz > 0')
        err = True

    if nbi['widy'] < 0:
        error('Invalid widy. Expected widy > 0')
        err = True

    if nbi['ashape'] not in [1, 2]:
        error('Invalid aperture shape. Expected 1 (rectangular) or 2 (circular)')
        err = True

    w = nbi['awidy'] < 0
    nw = len(nbi['awidy'][w])
    if nw > 0:
        error('Invalid awidy. Expected awidy >= 0.0')
        err = True

    w = nbi['awidz'] < 0
    nw = len(nbi['awidz'][w])
    if nw > 0:
        error('Invalid awidz. Expected awidz >= 0.0')
        err = True

    origin = inputs['origin']
    uvw_src = nbi['src']
    uvw_axis = nbi['axis']
    uvw_pos = uvw_src + nbi['adist'][0] * uvw_axis

    xyz_src = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_src, origin=origin)
    xyz_axis = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_axis)
    xyz_pos = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_pos, origin=origin)
    xyz_center = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], [0., 0., 0.], origin=origin)

    dis = np.sqrt(np.sum((xyz_src - xyz_pos) ** 2.))
    alpha = np.arctan2((xyz_pos[1] - xyz_src[1]), (xyz_pos[0] - xyz_src[0]))
    if dis != 0.:
        beta = np.arcsin((xyz_src[2] - xyz_pos[2]) / dis)
    else:
        beta = np.nan

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
    print('alpha = {} deg.'.format(alpha / np.pi * 180.))
    print('beta = {} deg.'.format(beta / np.pi * 180.))

    # Calculate grid center rc and sides length dr
    dr = np.array([inputs['xmax'] - inputs['xmin'], inputs['ymax'] - inputs['ymin'], inputs['zmax'] - inputs['zmin']], dtype=np.float64)
    rc = np.array([inputs['xmin'], inputs['ymin'], inputs['zmin']], dtype=np.float64) + 0.5 * dr

    # Check if beam centerline intersects beam grid
    length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_src, xyz_axis)

    print('Beam centerline - grid intersection length')
    print(length)
    if length <= 10.0:
        error('Beam centerline does not intersect grid')
        err = True

    if err:
        error('Invalid beam geometry. Exiting...', halt=True)
    else:
        success('Beam geometry is valid')

    return nbi
