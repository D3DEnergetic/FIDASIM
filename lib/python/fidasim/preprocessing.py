#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

#+#PREFIDA Source
#+ This file contains the source code for PREFIDA
#+***
from __future__ import print_function
import os
import numpy as np
import datetime
import h5py

from fidasim.utils import *

def check_dict_schema(schema, dic, desc=None):
    """
    #+#check_dict_schema
    #+ Check dict `dec` is formatted according to `schema`
    #+***
    #+##Input Arguments
    #+     **schema**: dict schema
    #+
    #+     **dic**: dict to check
    #+
    #+##Output Arguments
    #+     **err**: error code
    #+
    #+##Keyword Arguments
    #+     **desc**: description of dict `dic`
    #+
    #+##Example usage
    #+```python
    #+>>> dic = {'a':0, 'b':[1.d0,2.d0], 'c':"example"}
    #+>>> schema = {'a':{'dims':0,'type':[int]}, 'b':{'dims':[2],'type':[float, np.float64]}, 'c':{'dims':0,'type':[str]}  }
    #+
    #+>>> err = check_dict_schema(schema, dic, desc="Example dict")
    #+>>> print(err)
    #+    False
    #+```
    """
    if desc is None:
        desc = 'dict'

    err = False
    schema_keys = list(schema.keys())
    dic_keys = list(dic.keys())

    # Note extra variables
    for key in dic_keys:
        if key not in schema_keys:
            info('Extra variable "{}" found in "{}"'.format(key, desc))

    for key in schema_keys:
        # Note missing data
        if key not in dic_keys:
            error('"{}" is missing from "{}"'.format(key, desc))
            err = True
        else:
            # Check type
            if (schema[key]['dims'] == 0):
                if not isinstance(dic[key], tuple(schema[key]['type'])):
                    error('"{}" has the wrong type of {}. Expected {}'.format(key, type(dic[key]), schema[key]['type']))
                    err = True
            elif dic[key].dtype.type not in schema[key]['type']:
                error('"{}" has the wrong type of {}. Expected {}'.format(key, dic[key].dtype.type, schema[key]['type']))
                err = True

            # Check for NaNs or Inf
            if (not isinstance(dic[key], (str, dict, float, int))) and (np.bytes_ not in schema[key]['type']):
                if (dic[key][np.isnan(dic[key])].size > 0) or (dic[key][np.isinf(dic[key])].size > 0):
                    error('NaN or Infinity detected in "{}"'.format(key))
                    err = True

                # Check shape
                if not np.array_equal(dic[key].shape, schema[key]['dims']):
                    error('"{}" has the wrong shape of {}. Expected ({})'.format(key, dic[key].shape, schema[key]['dims']))
                    print('ndim({}) = {}'.format(key, dic[key].ndim))
                    err = True

            # Check shape
            if isinstance(dic[key], (str, int, float)):
                if (schema[key]['dims'] != 0) and (schema[key]['dims'] != [0]):
                    error('"{}" has the wrong shape. Expected ({})'.format(key, schema[key]['dims']))
                    err = True

    return err

def check_detector_overlap(d_cent, d_redge, d_tedge, d_shape, system_name, channel_ids):
    """
    Check if any detectors overlap within a diagnostic system.

    Input Arguments:
        d_cent: Detector center positions [3, nchan]
        d_redge: Detector radial edge positions [3, nchan]
        d_tedge: Detector toroidal edge positions [3, nchan]
        d_shape: Detector shape (1=rectangular, 2=circular) [nchan]
        system_name: Name of the diagnostic system (e.g., "NPA", "NC")
        channel_ids: Channel identifiers [nchan]

    Returns:
        overlaps: List of tuples (i, j) of overlapping detector pairs
    """
    nchan = d_cent.shape[1]
    overlaps = []

    # Check all pairs of detectors
    for i in range(nchan):
        for j in range(i + 1, nchan):
            # Calculate detector properties for both detectors
            cent_i = d_cent[:, i]
            cent_j = d_cent[:, j]

            # Calculate edge vectors
            e1_i = d_redge[:, i] - cent_i
            e2_i = d_tedge[:, i] - cent_i
            e1_j = d_redge[:, j] - cent_j
            e2_j = d_tedge[:, j] - cent_j

            # Calculate normal vectors
            normal_i = np.cross(e1_i, e2_i)
            normal_j = np.cross(e1_j, e2_j)

            # Normalize the normals
            normal_i_norm = normal_i / np.linalg.norm(normal_i)
            normal_j_norm = normal_j / np.linalg.norm(normal_j)

            # Check if detectors are roughly co-planar (parallel normals)
            dot_product = abs(np.dot(normal_i_norm, normal_j_norm))
            are_parallel = dot_product > 0.99  # Nearly parallel (within ~8 degrees)

            if are_parallel:
                # Check if they're in the same plane (distance from one center to the other's plane)
                displacement = cent_j - cent_i
                distance_to_plane = abs(np.dot(displacement, normal_i_norm))

                # Calculate effective radii (maximum distance from center to detector edge)
                if d_shape[i] == 1:  # Rectangular
                    radius_i = max(np.linalg.norm(e1_i), np.linalg.norm(e2_i))
                else:  # Circular
                    radius_i = np.linalg.norm(e1_i)

                if d_shape[j] == 1:  # Rectangular
                    radius_j = max(np.linalg.norm(e1_j), np.linalg.norm(e2_j))
                else:  # Circular
                    radius_j = np.linalg.norm(e1_j)

                # If in same plane, check if centers are close enough to overlap
                if distance_to_plane < 1.0:  # Within 1 cm of same plane
                    center_distance = np.linalg.norm(cent_j - cent_i)
                    # Conservative check: if distance between centers is less than sum of radii
                    if center_distance < (radius_i + radius_j):
                        overlaps.append((i, j))

    # Report overlaps if found
    if len(overlaps) > 0:
        warn('Overlapping detectors detected in {} diagnostic'.format(system_name))
        print('The following detector pairs may overlap:')
        for i, j in overlaps:
            id_i = channel_ids[i].decode('utf-8') if isinstance(channel_ids[i], bytes) else str(channel_ids[i])
            id_j = channel_ids[j].decode('utf-8') if isinstance(channel_ids[j], bytes) else str(channel_ids[j])
            print('  Channel {} (index {}) overlaps with Channel {} (index {})'.format(id_i, i, id_j, j))
        print('Please check your diagnostic geometry definition.')

    return overlaps

def check_chord_overlap(lens, axis, spot_size, system_name, channel_ids):
    """
    Check if any spectral chords overlap within a diagnostic system.

    Input Arguments:
        lens: Lens positions [3, nchan]
        axis: Chord axis directions [3, nchan]
        spot_size: Observation spot sizes [nchan]
        system_name: Name of the diagnostic system (e.g., "SPECTRAL", "FIDA")
        channel_ids: Channel identifiers [nchan]

    Returns:
        overlaps: List of tuples (i, j) of overlapping chord pairs
    """
    nchan = lens.shape[1]
    overlaps = []

    # Check all pairs of chords
    for i in range(nchan):
        for j in range(i + 1, nchan):
            # Calculate properties for both chords
            lens_i = lens[:, i]
            lens_j = lens[:, j]
            axis_i = axis[:, i]
            axis_j = axis[:, j]

            # Check if axes are nearly parallel
            dot_product = abs(np.dot(axis_i, axis_j))
            are_parallel = dot_product > 0.99  # Nearly parallel (within ~8 degrees)

            if are_parallel:
                # Check distance between lens positions
                lens_distance = np.linalg.norm(lens_j - lens_i)

                # Use spot sizes to determine potential overlap
                # If both spot sizes are defined, use their sum
                if spot_size[i] > 0 and spot_size[j] > 0:
                    threshold = spot_size[i] + spot_size[j]
                else:
                    # Use a default threshold of 2 cm if spot sizes not specified
                    threshold = 2.0

                if lens_distance < threshold:
                    overlaps.append((i, j))

    # Report overlaps if found
    if len(overlaps) > 0:
        warn('Overlapping chords detected in {} diagnostic'.format(system_name))
        print('The following chord pairs may overlap (same/similar lens position and parallel axes):')
        for i, j in overlaps:
            id_i = channel_ids[i].decode('utf-8') if isinstance(channel_ids[i], bytes) else str(channel_ids[i])
            id_j = channel_ids[j].decode('utf-8') if isinstance(channel_ids[j], bytes) else str(channel_ids[j])
            print('  Channel {} (index {}) overlaps with Channel {} (index {})'.format(id_i, i, id_j, j))
        print('Please check your diagnostic geometry definition.')

    return overlaps

def check_inputs(inputs, use_abs_path=True):
    """
    #+#check_inputs
    #+Checks if input dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: input dictionary
    #+
    #+##Output Arguments
    #+     **inputs**: Updated inputs dictionary
    #+
    #+##Example Usage
    #+```dist
    #+>>> inputs = check_inputs(inputs)
    #+```
    """
    info('Checking simulation settings...')
    err = False

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_int = {'dims': 0,
                'type': [int, np.int32, np.int64]}

    zero_long = {'dims': 0,
                 'type': [int, np.int32, np.int64]}

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
              'current_fractions': three_double,
              'pinj': zero_double,
              'einj': zero_double,
              'n_fida': zero_long,
              'n_nbi': zero_long,
              'n_pfida': zero_long,
              'n_pnpa': zero_long,
              'n_dcx': zero_long,
              'n_npa': zero_long,
              'n_halo': zero_long,
              'n_birth': zero_long,
              'ne_wght': zero_int,
              'np_wght': zero_int,
              'nphi_wght': zero_int,
              'emax_wght': zero_double,
              'ne_nc': zero_int,
              'np_nc': zero_int,
              'emax_nc_wght': zero_double,
              'nenergy_nc': zero_int,
              'emin_nc': zero_double,
              'emax_nc': zero_double,
              'nlambda_wght': zero_int,
              'lambdamin_wght': zero_double,
              'lambdamax_wght': zero_double,
              'calc_npa': zero_int,
              'calc_fida': zero_int,
              'calc_pnpa': zero_int,
              'calc_pfida': zero_int,
              'calc_bes': zero_int,
              'calc_dcx': zero_int,
              'calc_halo': zero_int,
              'calc_cold': zero_int,
              'calc_brems': zero_int,
              'calc_birth': zero_int,
              'calc_fida_wght': zero_int,
              'calc_npa_wght': zero_int,
              'calc_nc_wght': zero_int,              
              'calc_neutron': zero_int,
              'calc_neut_spec': zero_int,
              'calc_cfpd': zero_int,
              'calc_res': zero_int,
              'adaptive': zero_int,
              'split_tol': zero_double,
              'max_cell_splits': zero_int,
              'max_crossings': zero_int}

    # If user doesn't provide adaptive time step parameters, set default to off
    inputs.setdefault('adaptive', 0)
    inputs.setdefault('split_tol', 0.0)
    inputs.setdefault('max_cell_splits', 1)
    inputs.setdefault('max_crossings', 2)

    # Determine if this is a passive-only run
    is_passive_only = (inputs.get('calc_npa', 0) == 0 and
                       inputs.get('calc_fida', 0) == 0 and
                       inputs.get('calc_bes', 0) == 0 and
                       inputs.get('calc_dcx', 0) == 0 and
                       inputs.get('calc_halo', 0) == 0 and
                       inputs.get('calc_birth', 0) == 0 and
                       inputs.get('calc_nbi', 0) == 0 and
                       inputs.get('calc_fida_wght', 0) == 0 and
                       inputs.get('calc_npa_wght', 0) == 0)

    # For passive-only runs, beam parameters are optional
    if is_passive_only:
        # Create a modified schema without beam parameters
        passive_schema = schema.copy()
        # Remove beam-specific required parameters
        beam_params = ['nx', 'ny', 'nz', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax',
                      'alpha', 'beta', 'gamma', 'origin', 'current_fractions', 'pinj', 'einj']
        for param in beam_params:
            if param in passive_schema:
                del passive_schema[param]
        err = check_dict_schema(passive_schema, inputs, desc="simulation settings")
    else:
        err = check_dict_schema(schema, inputs, desc="simulation settings")

    if err:
        error('Invalid simulation settings. Exiting...', halt=True)

    # Normalize File Paths
    if use_abs_path:
        inputs['result_dir'] = os.path.abspath(inputs['result_dir'])

    # Only validate beam parameters if not in passive-only mode
    if not is_passive_only:
        if (inputs['alpha'] > 2. * np.pi) or (inputs['beta'] > 2. * np.pi) or (inputs['gamma'] > 2. * np.pi):
            error('Angles must be in radians')
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

    # These validations apply to both active and passive runs
    if inputs['lambdamin'] >= inputs['lambdamax']:
        error('Invalid wavelength range. Expected lambdamin < lamdbdamax')
        err = True

    if inputs['lambdamin_wght'] >= inputs['lambdamax_wght']:
        error('Invalid wavelength range. Expected lambdamin_wght < lamdbdamax_wght')
        err = True

    if inputs['emin_nc'] >= inputs['emax_nc']:
        error('Invalid neutron energy range. Expected emin_nc < emax_nc')
        err = True

    if (inputs['adaptive'] < 0):
        error('Invalid adaptive switch. Expected value >= 0')
        print('adaptive = {}'.format(inputs['adaptive']))
        err = True

    if (inputs['split_tol'] < 0.0) or (inputs['split_tol'] > 1.0):
        error('Invalid split tolerance. Expected value between 0 and 1')
        print('split_tol = {}'.format(inputs['split_tol']))
        err = True

    if (inputs['max_cell_splits'] < 1):
        error('Invalid max cell splits. Expected value >= 1')
        print('max_cell_splits = {}'.format(inputs['max_cell_splits']))
        err = True

    if (inputs['max_crossings'] < 2):
        error('Invalid max_crossings. Expected value >= 2')
        print('max_crossings = {}'.format(inputs['max_crossings']))
        err = True

    ps = os.path.sep
    input_file = inputs['result_dir'].rstrip(ps) + ps + inputs['runid'] + '_inputs.dat'
    equilibrium_file = inputs['result_dir'].rstrip(ps) + ps + inputs['runid'] + '_equilibrium.h5'
    geometry_file = inputs['result_dir'].rstrip(ps) + ps + inputs['runid'] + '_geometry.h5'
    distribution_file = inputs['result_dir'].rstrip(ps) + ps + inputs['runid'] + '_distribution.h5'
    neutrals_file = inputs['result_dir'].rstrip(ps) + ps + inputs['runid'] + '_neutrals.h5'

    inputs['input_file'] = input_file
    inputs['equilibrium_file'] = equilibrium_file
    inputs['geometry_file'] = geometry_file
    inputs['distribution_file'] = distribution_file
    inputs.setdefault('load_neutrals', 0)
    inputs.setdefault('output_neutral_reservoir', 1)
    inputs.setdefault('stark_components', 0)
    inputs.setdefault('flr', 2)
    inputs.setdefault('seed', -1)
    inputs.setdefault('verbose',1)
    inputs.setdefault('neutrals_file',neutrals_file)

    if err:
        error('Invalid simulation settings. Exiting...', halt=True)
    else:
        success('Simulation settings are valid')

    return inputs

def check_grid(grid):
    """
    #+#check_grid
    #+Checks if interpolation grid structure is valid
    #+***
    #+##Input Arguments
    #+     **grid**: Interpolation grid structure
    #+
    #+##Example Usage
    #+```python
    #+>>> check_grid(grid)
    #+```
    """
    err = False
    info('Checking interpolation grid...')

    if 'nr' not in grid:
        error('"nr" is missing from the interpolation grid')
        error('Invalid interpolation grid. Exiting...', halt=True)

    if 'nz' not in grid:
        error('"nz" is missing from the interpolation grid')
        error('Invalid interpolation grid. Exiting...', halt=True)

    if 'nphi' not in grid:
        info('"nphi" is missing from the interpolation grid, assuming axisymmetry')
        nphi = 1
    else:
        nphi = grid['nphi']
        # Check if phi array exists (for 3D geometry)
        if 'phi' in grid and nphi < 3:
            error('"phi" must have at least 3 elements when phi array is specified')
            error('Invalid interpolation grid. Exiting...', halt=True)
        elif 'phi' not in grid and nphi > 1:
            info('nphi={} specified without phi array - will be used for passive grid only'.format(nphi))


    nr = grid['nr']
    nz = grid['nz']

    zero_int = {'dims': 0,
                'type': [int, np.int32, np.int64]}

    nrnz_doub = {'dims': [nr, nz],
                 'type': [float, np.float64]}

    schema = {'nr': zero_int,
              'nz': zero_int,
              'r2d': nrnz_doub,
              'z2d': nrnz_doub,
              'r': {'dims': [nr],
                    'type': [float, np.float64]},
              'z': {'dims': [nz],
                    'type': [float, np.float64]}}

    if 'nphi' in grid:
        schema.setdefault('nphi',zero_int)
        # Only include phi in schema if it exists in the grid
        if 'phi' in grid:
            schema.setdefault('phi', {'dims': [nphi],
                                      'type': [float, np.float64]})

    err = check_dict_schema(schema, grid, desc="interpolation grid")
    if err:
        error('Invalid interpolation grid. Exiting...', halt=True)

    if not np.array_equal(grid['r'], np.sort(grid['r'])):
        error('r is not in ascending order')
        err = True

    if not np.array_equal(grid['z'], np.sort(grid['z'])):
        error('z is not in ascending order')
        err = True

    if not np.array_equal(grid['r'], grid['r2d'][:, 0]):
        error('r2d is defined incorrectly. Expected r == r2d[:, 0]')
        err = True

    if not np.array_equal(grid['z'], grid['z2d'][0, :]):
        error('z2d is defined incorrectly. Expected z == z2d[0, :]')
        err = True

    if err:
        error('Invalid interpolation grid. Exiting...', halt=True)
    else:
        success('Interpolation grid is valid')

        # Calculate and display grid resolution
        nr = grid['nr']
        nz = grid['nz']
        dr = (grid['r'][-1] - grid['r'][0]) / (nr - 1) if nr > 1 else 0.0
        dz = (grid['z'][-1] - grid['z'][0]) / (nz - 1) if nz > 1 else 0.0

        info('Grid resolution:')
        print('  Nr = {}, dR = {:.3f} cm'.format(nr, dr))
        print('  Nz = {}, dZ = {:.3f} cm'.format(nz, dz))

        if 'phi' in grid and 'nphi' in grid:
            nphi = grid['nphi']
            dphi = (grid['phi'][-1] - grid['phi'][0]) / (nphi - 1) if nphi > 1 else 0.0
            dphi_deg = np.degrees(dphi)
            print('  Nphi = {}, dPhi = {:.5f} rad ({:.2f} degrees)'.format(nphi, dphi, dphi_deg))
        elif 'nphi' in grid and grid['nphi'] > 1:
            info('  Nphi = {} specified without phi array - will be used for passive grid only'.format(grid['nphi']))

def check_beam(inputs, nbi):
    """
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
    """
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
        nbi['ashape'] = np.array(nbi['ashape'], dtype=int, ndmin=1)
        nbi['awidy'] = np.array(nbi['awidy'], dtype=float, ndmin=1)
        nbi['awidz'] = np.array(nbi['awidz'], dtype=float, ndmin=1)
        nbi['aoffy'] = np.array(nbi['aoffy'], dtype=float, ndmin=1)
        nbi['aoffz'] = np.array(nbi['aoffz'], dtype=float, ndmin=1)
        nbi['adist'] = np.array(nbi['adist'], dtype=float, ndmin=1)

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

    if nbi['naperture'] > 0:
        for aper in nbi['ashape']:
            if aper not in [1, 2]:
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

    # Machine coordinates
    origin = inputs['origin']
    uvw_src = nbi['src']
    uvw_axis = nbi['axis']
    if nbi['naperture'] == 0:
        uvw_pos = uvw_src + 100 * uvw_axis
    else:
        uvw_pos = uvw_src + nbi['adist'][0]

    # Convert to beam coordinates
    xyz_src = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_src, origin=origin)
    xyz_axis = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_axis)
    xyz_pos = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_pos, origin=origin)
    xyz_center = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], [0., 0., 0.], origin=origin)

    dis = np.sqrt(np.sum((xyz_src - xyz_pos) ** 2.))
    alpha = np.arctan2((xyz_pos[1] - xyz_src[1]), (xyz_pos[0] - xyz_src[0]))
    beta = np.arcsin((xyz_src[2] - xyz_pos[2]) / dis)

    print('Machine center in beam grid coordinates')
    print(xyz_center)
    print('Beam injection start point in machine coordinates')
    print(uvw_src)
    print('Beam injection start point in beam grid coordinates')
    print(xyz_src)

    if nbi['naperture'] != 0:
        print('First aperture position in machine coordinates')
        print(uvw_pos)
        print('First aperture position in beam grid coordinates')
        print(xyz_pos)
    else:
        print('Position of point 100cm along beam centerline in machine coordinates')
        print(uvw_pos)
        print('Position of point 100cm along beam centerline in beam grid coordinates')
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

def check_plasma(inputs, grid, plasma):
    """
    #+#check_plasma
    #+Checks if plasma paramters dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: Input dictionary
    #+
    #+     **grid**: Interpolation grid dictionary
    #+
    #+     **plasma**: Plasma parameters dictionary
    #+
    #+##Output Arguments
    #+     **plasma**: Updated plasma dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> plasma = check_plasma(inputs, grid, plasma)
    #+```
    """
    err = False
    info('Checking plasma parameters...')

    nr = grid['nr']
    nz = grid['nz']

    # Check for phi array to determine if we're truly 3D
    # nphi without phi array means passive grid only (plasma stays 2D)
    is_3d = 'phi' in grid and 'nphi' in grid and grid['nphi'] > 1

    if is_3d:
        nphi = grid['nphi']
    else:
        nphi = 1  # Treat as 2D for plasma arrays

    nthermal = plasma['nthermal']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_double = {'dims': 0,
                   'type': [float, np.float64]}

    zero_int = {'dims':0,
                'type': [int, np.int16, np.int32, np.int64]}

    nthermal_double = {'dims': [nthermal],
                       'type':[float, np.float64]}

    if is_3d:
        nrnznphi_double = {'dims': [nr, nz, nphi],
                           'type': [float, np.float64]}

        ntnrnznphi_double = {'dims': [nthermal, nr, nz, nphi],
                             'type': [float, np.float64]}

        nrnznphi_int = {'dims': [nr, nz, nphi],
                        'type': [np.int64]}
    else:
        nrnznphi_double = {'dims': [nr, nz],
                           'type': [float, np.float64]}

        ntnrnznphi_double = {'dims': [nthermal, nr, nz],
                             'type': [float, np.float64]}

        nrnznphi_int = {'dims': [nr, nz],
                        'type': [np.int64]}

    schema = {'time': zero_double,
              'nthermal': zero_int,
              'impurity_charge':zero_int,
              'species_mass': nthermal_double,
              'vr': nrnznphi_double,
              'vt': nrnznphi_double,
              'vz': nrnznphi_double,
              'dene': nrnznphi_double,
              'denimp': nrnznphi_double,
              'deni': ntnrnznphi_double,
              'denn': nrnznphi_double,
              'ti': nrnznphi_double,
              'te': nrnznphi_double,
              'zeff': nrnznphi_double,
              'mask': nrnznphi_int,
              'data_source': zero_string}

    err = check_dict_schema(schema, plasma, desc="plasma parameters")
    if err:
        error('Invalid plasma parameters. Exiting...', halt=True)

    if plasma['data_source'] == '':
        error('Invalid data source. An empty string is not a data source.')
        err = True

    # Electron density
    w = (plasma['dene'] < 0.)
    plasma['dene'][w] = 0.

    # Impurity density
    w = (plasma['denimp'] < 0.)
    plasma['denimp'][w] = 0.

    # Ion density
    w = (plasma['deni'] < 0.)
    plasma['deni'][w] = 0.

    # Zeff
    plasma['zeff'] = np.clip(plasma['zeff'],1.0,plasma['impurity_charge'])

    # Electron temperature
    w = (plasma['te'] < 0.)
    plasma['te'][w] = 0.

    # Ion temperature
    w = (plasma['ti'] < 0.)
    plasma['ti'][w] = 0.

    if (np.abs(plasma['time'] - inputs['time']) > 0.02):
        warn('Plasma time and input time do not match')
        print('Input time: ', inputs['time'])
        print('Plasma time: ', plasma['time'])

    # Add grid elements to plasma dict
    plasma.update(grid.copy())

    if err:
        error('Invalid plasma parameters. Exiting...', halt=True)
    else:
        success('Plasma parameters are valid')

    return plasma

def check_fields(inputs, grid, fields):
    """
    #+#check_fields
    #+Checks if electromagnetic fields dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: Input dictionary
    #+
    #+     **grid**: Interpolation grid dictionary
    #+
    #+     **fields**: Electromagnetic fields dictionary
    #+
    #+##Output Arguments
    #+     **fields**: Updated fields dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> fields = check_fields(inputs, grid, fields)
    #+```
    """
    err = False
    info('Checking electromagnetic fields...')

    nr = grid['nr']
    nz = grid['nz']
    # Check for phi array to determine if we're truly 3D
    # nphi without phi array means passive grid only (fields stay 2D)
    is_3d = 'phi' in grid and 'nphi' in grid and grid['nphi'] > 1

    if is_3d:
        nphi = grid['nphi']
    else:
        nphi = 1  # Treat as 2D for field arrays

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_double = {'dims': 0,
                   'type': [float, np.float64]}

    if is_3d:
        nrnznphi_double = {'dims': [nr, nz, nphi],
                           'type': [float, np.float64, np.float64]}

        nrnznphi_int = {'dims': [nr, nz, nphi],
                        'type': [int, np.int32, np.int64]}
    else:
        nrnznphi_double = {'dims': [nr, nz],
                           'type': [float, np.float64]}

        nrnznphi_int = {'dims': [nr, nz],
                        'type': [int, np.int32]}

    schema = {'time': zero_double,
              'br': nrnznphi_double,
              'bt': nrnznphi_double,
              'bz': nrnznphi_double,
              'er': nrnznphi_double,
              'et': nrnznphi_double,
              'ez': nrnznphi_double,
              'mask': nrnznphi_int,
              'data_source': zero_string}

    err = check_dict_schema(schema, fields, desc="electromagnetic fields")
    if err:
        error('Invalid electromagnetic fields. Exiting...', halt=True)

    if fields['data_source'] == '':
        error('Invalid data source. An empty string is not a data source.')
        err = True

    if np.abs(fields['time'] - inputs['time']) > 0.02:
        warn('Electromagnetic fields time and input time do not match')
        print('Input time: {}'.format(inputs['time']))
        print('Electromagnetic fields time: {}'.format(fields['time']))

    # Add grid elements to fields dict
    fields.update(grid.copy())

    if err:
        error('Invalid electromagnetic fields. Exiting...', halt=True)
    else:
        success('Electromagnetic fields are valid')

    return fields

def check_distribution(inputs, grid, dist):
    """
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
    """
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

        # Check for phi array to determine if we're truly 3D
        # nphi without phi array means passive grid only (distribution stays 2D)
        is_3d = 'phi' in grid and 'nphi' in grid and grid['nphi'] > 1

        if is_3d:
            nphi = grid['nphi']
        else:
            nphi = 1  # Treat as 2D for distribution arrays

        zero_string = {'dims': 0,
                       'type': [str]}

        zero_int = {'dims': 0,
                    'type': [int, np.int32, np.int64]}

        zero_double = {'dims': 0,
                       'type': [float, np.float64]}

        if is_3d:
            nrnznphi_double = {'dims': [nr, nz, nphi],
                               'type': [float, np.float64, np.float64]}
        else:
            nrnznphi_double = {'dims': [nr, nz],
                               'type': [float, np.float64]}

        schema = {'type': zero_int,
                  'nenergy': zero_int,
                  'npitch': zero_int,
                  'energy': {'dims': [nen],
                             'type': [float, np.float64]},
                  'pitch': {'dims': [npitch],
                            'type': [float, np.float64]},
                  'denf': nrnznphi_double,
                  'time': zero_double,
                  'data_source': zero_string}
        if is_3d:
            schema['f'] = {'dims': [nen, npitch, nr, nz, nphi],
                           'type': [float, np.float64]}
        else:
            schema['f'] = {'dims': [nen, npitch, nr, nz],
                           'type': [float, np.float64]}

        err = check_dict_schema(schema, dist, desc="fast-ion distribution")
        if err:
            error('Invalid fast-ion distribution. Exiting...', halt=True)

        # Add grid elements to plasma dict
        dist.update(grid.copy())

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

def check_spec(inputs, chords):
    """
    #+#check_spec
    #+Check if spectral geometry dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: input dictionary
    #+
    #+     **chords**: spectral geometry dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> check_spec(inputs, chords)
    #+```
    """
    err = False
    info('Checking FIDA/BES inputs...')

    chords_keys = list(chords.keys())

    if 'nchan' not in chords_keys:
        error('"nchan" is missing from the FIDA/BES geometry')
        err = True
        error('Invalid FIDA/BES geometry. Exiting...', halt=True)

    nchan = chords['nchan']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_long = {'dims': 0,
                 'type': [int, np.int32, np.int64]}

    nchan_double = {'dims': [nchan],
                    'type': [float, np.float64]}

    nchan_string = {'dims': [nchan],
                    'type': [np.bytes_, str]}

    three_nchan_float = {'dims': [3, nchan],
                         'type': [float, np.float64]}

    schema = {'data_source': zero_string,
              'nchan': zero_long,
              'system': zero_string,
              'id': nchan_string,
              'lens': three_nchan_float,
              'axis': three_nchan_float,
              'sigma_pi': nchan_double,
              'spot_size': nchan_double,
              'radius': nchan_double}

    err = check_dict_schema(schema, chords, desc="FIDA/BES geometry")
    if err:
        error('Invalid FIDA/BES geometry. Exiting...', halt=True)

    cross_arr = np.zeros(nchan, dtype=int)
    uvw_lens = chords['lens']
    uvw_axis = chords['axis']

    # Check if beam parameters are present (not passive-only mode)
    has_beam_params = all(key in inputs for key in ['alpha', 'beta', 'gamma', 'origin', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'])

    if has_beam_params:
        # ROTATE CHORDS INTO BEAM GRID COORDINATES
        xyz_lens = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_lens, origin=inputs['origin'])
        xyz_axis = uvw_to_xyz(inputs['alpha'], inputs['beta'], inputs['gamma'], uvw_axis)

        # Calculate grid center rc and sides length dr
        dr = np.array([inputs['xmax'] - inputs['xmin'], inputs['ymax'] - inputs['ymin'], inputs['zmax'] - inputs['zmin']], dtype=float)
        rc = np.array([inputs['xmin'], inputs['ymin'], inputs['zmin']], dtype=float) + 0.5 * dr

    for i in range(nchan):
        if not np.isclose(np.linalg.norm(uvw_axis[:, i]), 1, rtol=1e-03):
            error('Invalid optical axis for chord "' + str(chords['id'][i]) + '". Expected norm(axis) == 1')
            print(np.sum(uvw_axis[:, i] ** 2.))

        if has_beam_params:
            # Check if viewing chord intersects beam grid
            length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_lens[:, i], xyz_axis[:, i])
            if length <= 0.0:
                cross_arr[i] = 1

    if has_beam_params:
        wbad = (cross_arr == 1)
        nbad = cross_arr[wbad].size
        if nbad > 0:
            warn('The following {} chords do not cross the beam grid:'.format(nbad))
            warn('Chord ID: {}'.format(chords['id'][wbad]))

        wgood = (cross_arr == 0)
        ngood = cross_arr[wgood].size
        print('{} out of {} chords crossed the beam grid'.format(ngood, nchan))
        if ngood == 0:
            error('No channels intersect the beam grid')
            err = True
    else:
        # For passive-only runs, skip beam grid intersection check
        info('Passive-only mode: skipping beam grid intersection check for spectral chords')

    # Check for overlapping chords
    check_chord_overlap(chords['lens'], chords['axis'], chords['spot_size'], chords['system'], chords['id'])

    if err:
        error('Invalid FIDA/BES geometry. Exiting...', halt=True)
    else:
        success('FIDA/BES geometry is valid')

def check_npa(inp, npa):
    """
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
    """
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

    nchan_float = {'dims': [nchan],
                    'type': [float, np.float64]}

    nchan_string = {'dims': [nchan],
                    'type': [np.bytes_]}

    schema = {'data_source': zero_string,
              'nchan': zero_long,
              'system': zero_string,
              'id': nchan_string,
              'a_shape': nchan_int,
              'd_shape': nchan_int,
              'a_tedge': three_float,
              'a_redge': three_float,
              'a_cent': three_float,
              'd_tedge': three_float,
              'd_redge': three_float,
              'd_cent': three_float,
              'radius': nchan_float}

    err = check_dict_schema(schema, npa, desc="NPA geometry")
    if err:
        error('Invalid NPA geometry. Exiting...', halt=True)

    # Check detector/aperture shape
    w = np.logical_or(npa['d_shape'] > 2, npa['d_shape'] == 0)
    nw = len(npa['d_shape'][w])
    if nw != 0:
        error('Invalid detector shape. Expected 1 (rectagular) or 2 (circular)')
        print('Invalid indices: {}'.format(np.arange(len(npa['d_shape']))[w]))
        err = True

    w = np.logical_or(npa['a_shape'] > 2, npa['a_shape'] == 0)
    nw = len(npa['a_shape'][w])
    if nw != 0:
        error('Invalid aperture shape. Expected 1 (rectagular) or 2 (circular)')
        print('Invalid indices: {}'.format(np.arange(len(npa['a_shape']))[w]))
        err = True

    # Check if beam parameters are present (not passive-only mode)
    has_beam_params = all(key in inp for key in ['alpha', 'beta', 'gamma', 'origin', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'])

    err_arr = np.zeros(nchan, dtype=int)

    for i in range(nchan):
        uvw_det = npa['d_cent'][:, i]
        d_e1 = npa['d_redge'][:, i] - uvw_det
        d_e2 = npa['d_tedge'][:, i] - uvw_det

        uvw_aper = npa['a_cent'][:, i]
        a_e1 = npa['a_redge'][:, i] - uvw_aper
        a_e2 = npa['a_tedge'][:, i] - uvw_aper

        uvw_dir = uvw_aper - uvw_det

        if has_beam_params:
            # Calculate grid center rc and sides length dr (only once)
            if i == 0:
                dr = np.array([inp['xmax'] - inp['xmin'], inp['ymax'] - inp['ymin'], inp['zmax'] - inp['zmin']])
                rc = np.array([inp['xmin'], inp['ymin'], inp['zmin']]) + 0.5 * dr

            #Rotate chords into beam grid coordinates
            xyz_aper = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_aper, origin=inp['origin'])
            xyz_det = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_det, origin=inp['origin'])
            xyz_dir = xyz_aper - xyz_det
            xyz_dir = xyz_dir / np.sqrt(np.sum(xyz_dir * xyz_dir))

            # Check if npa chord intersects beam grid
            length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_det, xyz_dir)
            if length <= 0.:
                err_arr[i] = 1
            else:
                # Check if NPA detector is pointing in the right direction
                d_enter = np.sqrt(np.sum((r_enter - xyz_aper)**2))
                d_exit = np.sqrt(np.sum((r_exit - xyz_aper)**2))
                if d_exit < d_enter:
                    err_arr[i] = 1

        # Check that the detector and aperture point in the same direction
        d_e3 = np.cross(d_e1, d_e2)
        a_e3 = np.cross(a_e1, a_e2)
        a_dp = np.sum(uvw_dir*a_e3)
        d_dp = np.sum(uvw_dir*d_e3)
        dp = np.sum(d_e3 * a_e3)
        if (dp <= 0.) or (a_dp <= 0.) or (d_dp <= 0.):
            error('The detector and/or aperture plane normal vectors are pointing in the wrong direction. The NPA definition is incorrect.')
            err_arr[i] = 1

    if has_beam_params:
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
    else:
        # For passive-only mode, skip the beam grid intersection message
        info('Passive-only mode: skipping beam grid intersection check for NPA detectors')

    # Check for overlapping detectors
    check_detector_overlap(npa['d_cent'], npa['d_redge'], npa['d_tedge'], npa['d_shape'], 'NPA', npa['id'])

    if err:
        error('Invalid NPA geometry. Exiting...', halt=True)
    else:
        success('NPA geometry is valid')

def check_nc(inp, nc):
    """
    Checks if Neutron Collimator geometry structure is valid.

    Args:
        inp: Input dictionary.
        nc: Neutron Collimator geometry dictionary.

    Returns:
        None

    Raises:
        ValueError: If Neutron Collimator geometry is invalid, halts execution.
    """

    err_status = 0
    print("Checking Neutron Collimator geometry...")

    # Check if 'nchan' exists in nc
    if 'nchan' not in nc:
        raise ValueError('"nchan" is missing from the Neutron Collimator geometry')

    nchan = nc['nchan']

     # Define schema for dictionary structure validation, using similar setup as check_npa
    zero_string = {'dims': 0, 'type': [str]}
    zero_long = {'dims': 0, 'type': [int, np.int32, np.int64]}
    nchan_int = {'dims': [nchan], 'type': [int, np.int32, np.int64]}
    nchan_double = {'dims': [nchan], 'type': [float, np.float64]}
    nchan_string = {'dims': [nchan], 'type': [np.bytes_, str]}  # Allow both bytes and str
    three_nchan_float = {'dims': [3, nchan], 'type': [float, np.float64]}

    schema = {
        'data_source': zero_string,
        'nchan': zero_long,
        'system': zero_string,
        'id': nchan_string,
        'a_shape': nchan_int,  
        'd_shape': nchan_int,
        'a_tedge': three_nchan_float,
        'a_redge': three_nchan_float,
        'a_cent': three_nchan_float,
        'd_tedge': three_nchan_float,
        'd_redge': three_nchan_float,
        'd_cent': three_nchan_float,
        'radius': nchan_double
    }

    # Check if nc adheres to the schema
    err_status = check_dict_schema(schema, nc, desc="Neutron Collimator geometry")  # Assuming you have a check_dict_schema function
    if err_status == 1:
        raise ValueError('Invalid Neutron Collimator geometry. Exiting...')

    # Check detector/aperture shape
    invalid_detector_indices = np.where((nc['d_shape'] > 2) | (nc['d_shape'] == 0))[0]
    if len(invalid_detector_indices) > 0:
        raise ValueError(
            'Invalid detector shape. Expected 1 (rectangular) or 2 (circular)\n'
            f'Invalid indices: {invalid_detector_indices}'
        )

    invalid_aperture_indices = np.where((nc['a_shape'] > 2) | (nc['a_shape'] == 0))[0]
    if len(invalid_aperture_indices) > 0:
        raise ValueError(
            'Invalid aperture shape. Expected 1 (rectangular) or 2 (circular)\n'
            f'Invalid indices: {invalid_aperture_indices}'
        )
    """
    # Calculate grid center and side lengths
    dr = np.array([inp['xmax'] - inp['xmin'], inp['ymax'] - inp['ymin'], inp['zmax'] - inp['zmin']])
    rc = np.array([inp['xmin'], inp['ymin'], inp['zmin']]) + 0.5 * dr

    err_arr = np.zeros(nchan)
    err_arr = np.zeros(nchan)
    for i in range(nchan):
        # 1. Extract Detector and Aperture Information for Channel 'i'
        uvw_det = nc['d_cent'][:, i]  # Detector center in UVW coordinates
        d_e1 = nc['d_redge'][:, i] - uvw_det  # Detector edge vector 1
        d_e2 = nc['d_tedge'][:, i] - uvw_det  # Detector edge vector 2

        uvw_aper = nc['a_cent'][:, i]  # Aperture center in UVW coordinates
        a_e1 = nc['a_redge'][:, i] - uvw_aper  # Aperture edge vector 1
        a_e2 = nc['a_tedge'][:, i] - uvw_aper  # Aperture edge vector 2

        uvw_dir = uvw_aper - uvw_det  # Direction vector from detector to aperture

        # 2. Convert to XYZ Coordinates
        xyz_aper = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_aper, origin=inp['origin'])
        xyz_det = uvw_to_xyz(inp['alpha'], inp['beta'], inp['gamma'], uvw_det, origin=inp['origin'])
        xyz_dir = xyz_aper - xyz_det
        xyz_dir /= np.linalg.norm(xyz_dir)  # Normalize the direction vector

        # 3. Check for Intersection with Beam Grid
        length, r_enter, r_exit = aabb_intersect(rc, dr, xyz_det, xyz_dir)

        if length <= 0.0:
            err_arr[i] = 1  # Flag error if no intersection

        # 4. Check Detector Direction
        d_enter = np.linalg.norm(r_enter - xyz_aper)
        d_exit = np.linalg.norm(r_exit - xyz_aper)

        if d_exit < d_enter:
            err_arr[i] = 1  # Flag error if detector is pointing in the wrong direction

        # 5. Check Detector and Aperture Alignment
        d_e3 = np.cross(d_e1, d_e2)  # Detector normal vector
        a_e3 = np.cross(a_e1, a_e2)  # Aperture normal vector

        d_dp = np.dot(uvw_dir, d_e3)
        a_dp = np.dot(uvw_dir, a_e3)
        dp = np.dot(d_e3, a_e3)

        if (a_dp <= 0.0) or (d_dp <= 0.0) or (dp <= 0.0):
            raise ValueError(
                'The detector and/or aperture plane normal vectors are '
                'pointing in the wrong direction. The Neutron Collimator '
                'definition is incorrect.'
            )
            err_arr[i] = 1

        # Check how many channels crossed the beam grid
    crossed_grid_indices = np.where(err_arr == 0)[0]
    num_crossed = len(crossed_grid_indices)
    num_missed = len(np.where(err_arr != 0)[0])

    print(f"{num_crossed} out of {nchan} channels crossed the beam grid")

    if num_crossed == 0:
        raise ValueError('No channels intersect the beam grid')

    if num_missed > 0:
        print('Some channels did not intersect the beam grid')
        print(f'Number missed: {num_missed}')
        print('Missed channels:')
        print(f'    {nc["id"][np.where(err_arr != 0)[0]]}')  # Assuming 'id' is a numpy array
    """

    # Check for overlapping detectors
    check_detector_overlap(nc['d_cent'], nc['d_redge'], nc['d_tedge'], nc['d_shape'], 'NC', nc['id'])

    print('Neutron Collimator geometry is valid')


def check_cfpd(inp, cfpd):
    """
    #+#check_cfpd
    #+Checks if CFPD geometry dictionary is valid
    #+***
    #+##Input Arguments
    #+     **inputs**: input dictionary
    #+
    #+     **cfpd**: CFPD geometry dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> check_cfpd(inputs, cfpd)
    #+```
    """
    err = False
    info('Checking CFPD geometry...')

    cfpd_keys = cfpd.keys()

    if 'nchan' not in cfpd_keys:
        error('"nchan" is missing from the CFPD geometry')
        err = True
        error('Invalid CFPD geometry. Exiting...', halt=True)

    nchan = cfpd['nchan']

    zero_string = {'dims': 0,
                   'type': [str]}

    zero_long = {'dims': 0,
                 'type': [int, np.int32, np.int64]}

    three_float = {'dims': [3, nchan],
                   'type': [float, np.float64]}

    nchan_int = {'dims': [nchan],
                 'type': [int, np.int32, np.int64]}

    nchan_float = {'dims': [nchan],
                    'type': [float, np.float64]}

    nchan_string = {'dims': [nchan],
                    'type': [np.bytes_]}

    schema = {'data_source': zero_string,
              'nchan': zero_long,
              'system': zero_string,
              'id': nchan_string,
              'a_shape': nchan_int,
              'd_shape': nchan_int,
              'a_tedge': three_float,
              'a_redge': three_float,
              'a_cent': three_float,
              'd_tedge': three_float,
              'd_redge': three_float,
              'd_cent': three_float,
              'radius': nchan_float}

    err = check_dict_schema(schema, cfpd, desc="CFPD geometry")
    if err:
        error('Invalid CFPD geometry. Exiting...', halt=True)

    # Check detector/aperture shape
    w = np.logical_or(cfpd['d_shape'] > 2, cfpd['d_shape'] == 0)
    nw = len(cfpd['d_shape'][w])
    if nw != 0:
        error('Invalid detector shape. Expected 1 (rectagular) or 2 (circular)')
        print('Invalid indices: {}'.format(np.arange(len(cfpd['d_shape']))[w]))
        err = True

    w = np.logical_or(cfpd['a_shape'] > 2, cfpd['a_shape'] == 0)
    nw = len(cfpd['a_shape'][w])
    if nw != 0:
        error('Invalid aperture shape. Expected 1 (rectagular) or 2 (circular)')
        print('Invalid indices: {}'.format(np.arange(len(cfpd['a_shape']))[w]))
        err = True

    err_arr = np.zeros(nchan, dtype=int)
    for i in range(nchan):
        uvw_det = cfpd['d_cent'][:, i]
        d_e1 = cfpd['d_redge'][:, i] - uvw_det
        d_e2 = cfpd['d_tedge'][:, i] - uvw_det

        uvw_aper = cfpd['a_cent'][:, i]
        a_e1 = cfpd['a_redge'][:, i] - uvw_aper
        a_e2 = cfpd['a_tedge'][:, i] - uvw_aper

        uvw_dir = uvw_aper - uvw_det

        # Check that the detector and aperture point in the same direction
        d_e3 = np.cross(d_e1, d_e2)
        a_e3 = np.cross(a_e1, a_e2)
        a_dp = np.sum(uvw_dir*a_e3)
        d_dp = np.sum(uvw_dir*d_e3)
        dp = np.sum(d_e3 * a_e3)
        if (dp <= 0.) or (a_dp <= 0.) or (d_dp <= 0.):
            error('The detector and/or aperture plane normal vectors are pointing in the wrong direction. The CFPD definition is incorrect.')
            err_arr[i] = 1

    w = (err_arr == 0)
    nw = err_arr[w].size
    if nw == 0:
        err = True

    # Check for overlapping detectors
    check_detector_overlap(cfpd['d_cent'], cfpd['d_redge'], cfpd['d_tedge'], cfpd['d_shape'], 'CFPD', cfpd['id'])

    if err:
        error('Invalid CFPD geometry. Exiting...', halt=True)
    else:
        success('CFPD geometry is valid')

def write_namelist(filename, inputs):
    """
    #+#write_namelist
    #+Writes namelist file
    #+***
    #+##Input Arguments
    #+     **filename**: Name of the namelist file
    #+
    #+     **inputs**: Input dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> write_namelist(filename, inputs)
    #+```
    """
    info("Writing namelist file...")

    fidasim_version = get_version(get_fidasim_dir())

    with open(filename, "w") as f:
        f.write("!! Created: {}\n".format(datetime.datetime.now()))
        f.write("!! FIDASIM version: {}\n".format(fidasim_version))
        f.write("!! Comment: {}\n".format(inputs['comment']))
        f.write("&fidasim_inputs\n\n")

        f.write("!! Shot Info\n")
        f.write("shot = {:d}    !! Shot Number\n".format(inputs['shot']))
        f.write("time = {:f}    !! Time [s]\n".format(inputs['time']))
        f.write("runid = '{}'   !! runID\n".format(inputs['runid']))
        f.write("result_dir = '{}'    !! Result Directory\n\n".format(inputs['result_dir']))

        f.write("!! Input Files\n")
        f.write("tables_file = '{}'   !! Atomic Tables File\n".format(inputs['tables_file']))
        f.write("equilibrium_file = '" + inputs['equilibrium_file'] + "'    !! File containing plasma parameters and fields\n")
        f.write("geometry_file = '" + inputs['geometry_file'] + "'    !! File containing NBI and diagnostic geometry\n")
        f.write("distribution_file = '" + inputs['distribution_file'] + "'    !! File containing fast-ion distribution\n\n")

        f.write("!! Simulation Switches\n")
        f.write("calc_bes = {:d}    !! Calculate NBI Spectra\n".format(inputs['calc_bes']))
        f.write("calc_dcx = {:d}    !! Calculate Direct CX Spectra\n".format(inputs['calc_dcx']))
        f.write("calc_halo = {:d}    !! Calculate Halo Spectra\n".format(inputs['calc_halo']))
        f.write("calc_cold = {:d}    !! Calculate Cold D-alpha Spectra\n".format(inputs['calc_cold']))
        f.write("calc_brems = {:d}    !! Calculate Bremsstrahlung\n".format(inputs['calc_brems']))
        f.write("calc_fida = {:d}    !! Calculate FIDA Spectra\n".format(inputs['calc_fida']))
        f.write("calc_npa = {:d}   !! Calculate NPA\n".format(inputs['calc_npa']))
        f.write("calc_pfida = {:d}    !! Calculate Passive FIDA Spectra\n".format(inputs['calc_pfida']))
        f.write("calc_pnpa = {:d}   !! Calculate Passive NPA\n".format(inputs['calc_pnpa']))
        f.write("calc_neutron = {:d}   !! Calculate B-T Neutron Rate\n".format(inputs['calc_neutron']))
        f.write("calc_neut_spec = {:d}   !! Calculate Neutron Spectra\n".format(inputs['calc_neut_spec']))
        f.write("calc_cfpd = {:d}   !! Calculate B-T CFPD Energy Resolved Count Rate\n".format(inputs['calc_cfpd']))
        f.write("calc_birth = {:d}    !! Calculate Birth Profile\n".format(inputs['calc_birth']))
        f.write("calc_fida_wght = {:d}    !! Calculate FIDA weights\n".format(inputs['calc_fida_wght']))
        f.write("calc_npa_wght = {:d}    !! Calculate NPA weights\n".format(inputs['calc_npa_wght']))
        f.write("calc_nc_wght = {:d}    !! Calculate NC weights\n".format(inputs['calc_nc_wght']))
        f.write("calc_res = {:d}    !! Calculate spatial resolution\n".format(inputs['calc_res']))

        f.write("\n!! Advanced Settings\n")
        f.write("seed = {:d}    !! RNG Seed. If seed is negative a random seed is used\n".format(inputs['seed']))
        f.write("flr = {:d}    !! Turn on Finite Larmor Radius corrections\n".format(inputs['flr']))
        f.write("load_neutrals = {:d}    !! Load neutrals from neutrals file\n".format(inputs['load_neutrals']))
        f.write("output_neutral_reservoir = {:d}    !! Output neutral reservoir to neutrals file\n".format(inputs['output_neutral_reservoir']))
        f.write("neutrals_file = '" + inputs['neutrals_file'] + "'    !! File containing the neutral density\n")
        f.write("stark_components = {:d}    !! Output stark components\n".format(inputs['stark_components']))
        f.write("verbose = {:d}    !! Verbose\n\n".format(inputs['verbose']))

        f.write("!! Monte Carlo Settings\n")
        f.write("n_fida = {:d}    !! Number of FIDA mc particles\n".format(inputs['n_fida']))
        f.write("n_npa = {:d}    !! Number of NPA mc particles\n".format(inputs['n_npa']))
        f.write("n_pfida = {:d}    !! Number of Passive FIDA mc particles\n".format(inputs['n_pfida']))
        f.write("n_pnpa = {:d}    !! Number of Passive NPA mc particles\n".format(inputs['n_pnpa']))
        f.write("n_nbi = {:d}    !! Number of NBI mc particles\n".format(inputs['n_nbi']))
        f.write("n_halo = {:d}    !! Number of HALO mc particles\n".format(inputs['n_halo']))
        f.write("n_dcx = {:d}     !! Number of DCX mc particles\n".format(inputs['n_dcx']))
        f.write("n_birth = {:d}    !! Number of BIRTH mc particles\n\n".format(inputs['n_birth']))

        # Check if this is a passive-only run
        is_passive_only = ('pinj' not in inputs or 'einj' not in inputs or
                          'current_fractions' not in inputs or 'nx' not in inputs)

        if not is_passive_only:
            f.write("!! Neutral Beam Settings\n")
            f.write("ab = {:f}     !! Beam Species mass [amu]\n".format(inputs['ab']))
            f.write("pinj = {:f}     !! Beam Power [MW]\n".format(inputs['pinj']))
            f.write("einj = {:f}     !! Beam Energy [keV]\n".format(inputs['einj']))
            f.write("current_fractions(1) = {:f} !! Current Fractions (Full component)\n".format(inputs['current_fractions'][0]))
            f.write("current_fractions(2) = {:f} !! Current Fractions (Half component)\n".format(inputs['current_fractions'][1]))
            f.write("current_fractions(3) = {:f} !! Current Fractions (Third component)\n\n".format(inputs['current_fractions'][2]))

            f.write("!! Beam Grid Settings\n")
            f.write("nx = {:d}    !! Number of cells in X direction (Into Plasma)\n".format(inputs['nx']))
            f.write("ny = {:d}    !! Number of cells in Y direction\n".format(inputs['ny']))
            f.write("nz = {:d}    !! Number of cells in Z direction\n".format(inputs['nz']))
            f.write("xmin = {:f}     !! Minimum X value [cm]\n".format(inputs['xmin']))
            f.write("xmax = {:f}     !! Maximum X value [cm]\n".format(inputs['xmax']))
            f.write("ymin = {:f}     !! Minimum Y value [cm]\n".format(inputs['ymin']))
            f.write("ymax = {:f}     !! Maximum Y value [cm]\n".format(inputs['ymax']))
            f.write("zmin = {:f}     !! Minimum Z value [cm]\n".format(inputs['zmin']))
            f.write("zmax = {:f}     !! Maximum Z value [cm]\n\n".format(inputs['zmax']))

            f.write("!! Tait-Bryan Angles for z-y`-x`` rotation\n")
            f.write("alpha = {:f}     !! Rotation about z-axis [rad]\n".format(inputs['alpha']))
            f.write("beta  = {:f}     !! Rotation about y`-axis [rad]\n".format(inputs['beta']))
            f.write("gamma = {:f}     !! Rotation about x``-axis [rad]\n\n".format(inputs['gamma']))

            f.write("!! Beam Grid origin in machine coordinates (cartesian)\n")
            f.write("origin(1) = {:f}     !! U value [cm]\n".format(inputs['origin'][0]))
            f.write("origin(2) = {:f}     !! V value [cm]\n".format(inputs['origin'][1]))
            f.write("origin(3) = {:f}     !! W value [cm]\n\n".format(inputs['origin'][2]))
        else:
            f.write("!! Passive-only run - no beam parameters needed\n")
            f.write("!! Fast-ion Species mass\n")
            f.write("ab = {:f}     !! Fast-ion Species mass [amu]\n\n".format(inputs['ab']))

        f.write("!! Wavelength Grid Settings\n")
        f.write("nlambda = {:d}    !! Number of Wavelengths\n".format(inputs['nlambda']))
        f.write("lambdamin = {:f}    !! Minimum Wavelength [nm]\n".format(inputs['lambdamin']))
        f.write("lambdamax = {:f}    !! Maximum Wavelength [nm]\n\n".format(inputs['lambdamax']))

        f.write("!! Weight Function Settings\n")
        f.write("ne_wght = {:d}    !! Number of Energies for Weights\n".format(inputs['ne_wght']))
        f.write("np_wght = {:d}    !! Number of Pitches for Weights\n".format(inputs['np_wght']))
        f.write("nphi_wght = {:d}    !! Number of Gyro-angles for Weights\n".format(inputs['nphi_wght']))
        f.write("emax_wght = {:f}    !! Maximum Energy for Weights [keV]\n".format(inputs['emax_wght']))
        f.write("nlambda_wght = {:d}    !! Number of Wavelengths for Weights \n".format(inputs['nlambda_wght']))
        f.write("lambdamin_wght = {:f}    !! Minimum Wavelength for Weights [nm]\n".format(inputs['lambdamin_wght']))
        f.write("lambdamax_wght = {:f}    !! Maximum Wavelength for Weights [nm]\n".format(inputs['lambdamax_wght']))
        f.write("ne_nc = {:d}    !! Number of Energies for NC Weights\n".format(inputs['ne_nc']))
        f.write("np_nc = {:d}    !! Number of Pitches for NC Weights\n".format(inputs['np_nc']))
        f.write("emax_nc_wght = {:f}    !! Maximum Energy for NC Weights [keV]\n\n".format(inputs['emax_nc_wght']))

        f.write("!! Neutron Spectra Settings\n")
        f.write("nenergy_nc = {:d}    !! Number of Energy Bins for Neutron Spectra\n".format(inputs['nenergy_nc']))
        f.write("emin_nc = {:f}    !! Minimum Neutron Energy [keV]\n".format(inputs['emin_nc']))
        f.write("emax_nc = {:f}    !! Maximum Neutron Energy [keV]\n\n".format(inputs['emax_nc']))

        f.write("!! Adaptive Time Step Settings\n")
        f.write("adaptive = {:d}    !! Adaptive switch, 0:split off, 1:dene, 2:denn, 3:denf, 4:deni, 5:denimp, 6:te, 7:ti\n".format(inputs['adaptive']))
        f.write("split_tol = {:f}    !! Tolerance for change in plasma parameter, number of cell splits is proportional to 1/split_tol\n".format(inputs['split_tol']))
        f.write("max_cell_splits = {:d}    !! Maximum number of times a cell can be split\n\n".format(inputs['max_cell_splits']))
        f.write("max_crossings = {:d}    !! Maximum number of times a neutral/LOS can cross the plasma boundary\n\n".format(inputs['max_crossings']))
        f.write("/\n\n")

    success("Namelist file created: {}\n".format(filename))

def write_geometry(filename, nbi=None, spec=None, npa=None, nc=None, cfpd=None):
    """
    #+#write_geometry
    #+Write geometry values to a HDF5 file
    #+***
    #+##Input Arguments
    #+     **filename**: Name of the geometry file
    #+
    #+##Keyword Arguments
    #+     **nbi**: Optional, NBI geometry structure (not needed for passive-only runs)
    #+
    #+     **spec**: Optional, Spectral geometry structure
    #+
    #+     **npa**: Optional, NPA geometry structure
    #+
    #+     **cfpd**: Optional, CFPD geometry structure
    #+
    #+##Example Usage
    #+```python
    #+>>> # With beam
    #+>>> write_geometry(filename, nbi=nbi, spec=spec, npa=npa, cfpd=cfpd)
    #+>>> # Passive-only (no beam)
    #+>>> write_geometry(filename, spec=spec, npa=npa, nc=nc)
    #+```
    """
    info('Writing geometry file...')

    # Create and open h5 file
    with h5py.File(filename, 'w') as hf:
        # File attributes
        hf.attrs['description'] = 'Geometric quantities for FIDASIM'

        # Create nbi group only if nbi geometry is provided
        if nbi is not None:
            g_nbi = hf.create_group('nbi')

            # nbi att
            g_nbi.attrs['description'] = 'Neutral Beam Geometry'
            g_nbi.attrs['coordinate_system'] = 'Right-handed cartesian'

            nbi_description = {'data_source': 'Source of the NBI geometry',
                               'name': 'Beam name',
                               'src': 'Position of the center of the beam source grid',
                               'axis':'Axis of the beam centerline: Centerline(t) = src + axis*t ',
                               'focy': 'Horizonal focal length of the beam',
                               'focz': 'Vertical focal length of the beam',
                               'divy': 'Horizonal divergences of the beam. One for each energy component',
                               'divz': 'Vertical divergences of the beam. One for each energy component',
                               'widy': 'Half width of the beam source grid',
                               'widz': 'Half height of the beam source grid',
                               'shape':'Shape of the beam source grid: 1="rectangular", 2="circular"',
                               'naperture': 'Number of apertures',
                               'ashape': 'Shape of the aperture(s): 1="rectangular", 2="circular"',
                               'awidy': 'Half width of the aperture(s)',
                               'awidz': 'Half height of the aperture(s)',
                               'aoffy': 'Horizontal (y) offset of the aperture(s) relative to the +x aligned beam centerline',
                               'aoffz': 'Vertical (z) offset of the aperture(s) relative to the +x aligned beam centerline',
                               'adist': 'Distance from the center of the beam source grid to the aperture(s) plane'}

            nbi_units = {'src': 'cm',
                         'axis': 'cm',
                         'focy': 'cm',
                         'focz': 'cm',
                         'divy': 'radians',
                         'divz': 'radians',
                         'widy': 'cm',
                         'widz': 'cm',
                         'awidy': 'cm',
                         'awidz': 'cm',
                         'aoffy': 'cm',
                         'aoffz': 'cm',
                         'adist': 'cm'}
            
            # Define dimension names for xarray/pandas compatibility
            # Note: Removing dimension labels for small configuration arrays to avoid xarray issues
            # These are parameter arrays, not data on coordinate grids             
            nbi_dim_names = {}

            write_data(g_nbi, nbi, desc = nbi_description, units=nbi_units,
                  dim_names=nbi_dim_names, name='nbi')

        if spec is not None:
            # Create spec group
            g_spec = hf.create_group('spec')

            # Spectroscopic attributes
            g_spec.attrs['description'] = 'FIDA/BES Chord Geometry'
            g_spec.attrs['coordinate_system'] = 'Right-handed cartesian'

            # Define description attributes
            spec_description = {'data_source': 'Source of the chord geometry',
                                'nchan': 'Number of channels',
                                'system': 'Names of the different spectrocopic systems',
                                'id': 'Line of sight ID',
                                'lens': 'Positions of the lenses',
                                'axis': 'Optical axis of the lines of sight: LOS(t) = lens + axis*t ',
                                'radius': 'Line of sight radius at midplane or tangency point',
                                'sigma_pi': 'Ratio of the intensities of the sigma and pi stark lines. Measured quantity',
                                'spot_size': 'Radius of spot size'}

            spec_units = {'lens': 'cm',
                          'axis': 'cm',
                          'radius': 'cm',
                          'spot_size': 'cm'}

            # Define dimension names for xarray/pandas compatibility
            # Note: Removing dimension labels for configuration arrays to avoid xarray issues
            spec_dim_names = {}

            write_data(g_spec, spec, desc=spec_description, units=spec_units,
                      dim_names=spec_dim_names, name='spec')

        if npa is not None:
            # Create npa group
            g_npa = hf.create_group('npa')

            # Group attributes
            g_npa.attrs['description'] = 'NPA Geometry'
            g_npa.attrs['coordinate_system'] = 'Right-handed cartesian'

            # Dataset attributes
            npa_description = {'data_source': 'Source of the NPA geometry',
                               'nchan': 'Number of channels',
                               'system': 'Names of the different NPA systems',
                               'id': 'Line of sight ID',
                               'd_shape': 'Shape of the detector: 1="rectangular", 2="circular"',
                               'd_cent': 'Center of the detector',
                               'd_tedge': 'Center of the detectors top edge',
                               'd_redge': 'Center of the detectors right edge',
                               'a_shape': 'Shape of the aperture: 1="rectangular", 2="circular"',
                               'a_cent': 'Center of the aperture',
                               'a_tedge': 'Center of the apertures top edge',
                               'a_redge': 'Center of the apertures right edge',
                               'radius': 'Line of sight radius at midplane or tangency point'}

            npa_units = {'d_cent': 'cm',
                         'd_tedge': 'cm',
                         'd_redge': 'cm',
                         'a_cent': 'cm',
                         'a_tedge': 'cm',
                         'radius': 'cm',
                         'a_redge': 'cm'}

            # Define dimension names for xarray/pandas compatibility
            # Note: Removing dimension labels for configuration arrays to avoid xarray issues
            npa_dim_names = {}

            write_data(g_npa, npa, desc = npa_description, units = npa_units,
                      dim_names=npa_dim_names, name='npa')
        
        if nc is not None:
            # Create npa group
            g_nc = hf.create_group('nc')

            # Group attributes
            g_nc.attrs['description'] = 'NC Geometry'
            g_nc.attrs['coordinate_system'] = 'Right-handed cartesian'

            # Dataset attributes
            nc_description = {'data_source': 'Source of the NC geometry',
                               'nchan': 'Number of channels',
                               'system': 'Names of the different NC systems',
                               'id': 'Line of sight ID',
                               'd_shape': 'Shape of the detector: 1="rectangular", 2="circular"',
                               'd_cent': 'Center of the detector',
                               'd_tedge': 'Center of the detectors top edge',
                               'd_redge': 'Center of the detectors right edge',
                               'a_shape': 'Shape of the aperture: 1="rectangular", 2="circular"',
                               'a_cent': 'Center of the aperture',
                               'a_tedge': 'Center of the apertures top edge',
                               'a_redge': 'Center of the apertures right edge',
                               'radius': 'Line of sight radius at midplane or tangency point'}

            nc_units = {'d_cent': 'cm',
                         'd_tedge': 'cm',
                         'd_redge': 'cm',
                         'a_cent': 'cm',
                         'a_tedge': 'cm',
                         'radius': 'cm',
                         'a_redge': 'cm'}

            # Define dimension names for xarray/pandas compatibility
            # Note: Removing dimension labels for configuration arrays to avoid xarray issues
            nc_dim_names = {}

            write_data(g_nc, nc, desc = nc_description, units = nc_units,
                      dim_names=nc_dim_names, name='nc')

        if cfpd is not None:
            # Create cfpd group
            g_cfpd = hf.create_group('cfpd')

            # Group attributes
            g_cfpd.attrs['description'] = 'CFPD Geometry'
            g_cfpd.attrs['coordinate_system'] = 'Right-handed cartesian'

            # Dataset attributes
            cfpd_description = {'data_source': 'Source of the CFPD geometry',
                               'nchan': 'Number of channels',
                               'nrays': 'Number of rays',
                               'nsteps': 'Maximum number of orbit steps',
                               'nenergy': 'Number of energies',
                               'nactual': 'Number of orbital spatial steps',
                               'system': 'Names of the different CFPD systems',
                               'id': 'Line of sight ID',
                               'd_shape': 'Shape of the detector: 1="rectangular", 2="circular"',
                               'd_cent': 'Center of the detector',
                               'd_tedge': 'Center of the detectors top edge',
                               'd_redge': 'Center of the detectors right edge',
                               'a_shape': 'Shape of the aperture: 1="rectangular", 2="circular"',
                               'a_cent': 'Center of the aperture',
                               'a_tedge': 'Center of the apertures top edge',
                               'a_redge': 'Center of the apertures right edge',
                               'radius': 'Line of sight radius at midplane or tangency point',
                               'earray': 'Energy array',
                               'sightline': 'Velocity and position in (R,Phi,Z)',
                               'daomega': 'Transmission factor'}

            cfpd_units = {'d_cent': 'cm',
                         'd_tedge': 'cm',
                         'd_redge': 'cm',
                         'a_cent': 'cm',
                         'a_tedge': 'cm',
                         'radius': 'cm',
                         'a_redge': 'cm',
                         'earray': 'keV',
                         'sightline': 'cm/s and cm',
                         'daomega': 'cm^2'}

            # Define dimension names for xarray/pandas compatibility
            # Note: Only keeping dimension labels where we have proper scales
            cfpd_dim_names = {'earray': ['energy']}

            cfpd_coord_scales = {'earray': 'energy'}

            write_data(g_cfpd, cfpd, desc = cfpd_description, units = cfpd_units,
                      dim_names=cfpd_dim_names, coord_scales=cfpd_coord_scales, name='cfpd')

    if os.path.isfile(filename):
        success('Geometry file created: ' + filename)
    else:
        error('Geometry file creation failed.')

def write_equilibrium(filename, plasma, fields):
    """
    #+#write_equilibrium
    #+Write MHD equilibrium values to a HDF5 file
    #+***
    #+##Input Arguments
    #+     **filename**: Name of the equilibrium file
    #+
    #+     **plasma**: Plasma dictionary
    #+
    #+     **fields**: Electromagnetic fields dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> write_equilibrium(filename, plasma, fields)
    #+```
    """
    info('Writing equilibrium file...')

    with h5py.File(filename, 'w') as hf:
        # File attribute
        hf.attrs['description'] = 'Plasma Parameters and Electromagnetic Fields for FIDASIM'

        # Create plasma group
        g_plasma = hf.create_group('plasma')

        # Plasma Attributes
        g_plasma.attrs['description'] = 'Plasma Parameters'
        g_plasma.attrs['coordinate_system'] = 'Cylindrical'

        # Dataset attributes
        plasma_description = {'data_source': 'Source of the plasma parameters',
                              'time': 'Time',
                              'nthermal': 'Number of hydrogenic, thermal ion species',
                              'dene': 'Electron Number Density: Dene(r,z,[phi])',
                              'denimp': 'Impurity Number Density: Denimp(r,z,[phi])',
                              'deni': 'Ion Number Density: Deni(species,r,z,[phi])',
                              'species_mass': 'Mass of ion species: species_mass(species)',
                              'impurity_charge': 'Charge number of main impurity species',
                              'te': 'Electron Temperature: Te(r,z,[phi])',
                              'ti': 'Ion Temperature: Ti(r,z,[phi])',
                              'zeff': 'Effective Nuclear Charge: Zeff(r,z,[phi])',
                              'denn': 'Cold/Edge neutral density: Denn(r,z,[phi])',
                              'vr': 'Bulk plasma flow in the r-direction: Vr(r,z,[phi])',
                              'vt': 'Bulk plasma flow in the theta/torodial-direction: Vt(r,z,[phi])',
                              'vz': 'Bulk plasma flow in the z-direction: Vz(r,z,[phi])',
                              'nr': 'Number of R values',
                              'nz': 'Number of Z values',
                              'r': 'Radius',
                              'z': 'Z',
                              'r2d': 'Radius grid: R(r,z)',
                              'z2d': 'Z grid: Z(r,z)',
                              'mask': 'Boolean mask that indicates where the plasma parameters are well defined'}

        plasma_units = {'time': 's',
                        'species_mass': 'amu',
                        'dene': 'cm^-3',
                        'deni': 'cm^-3',
                        'denn': 'cm^-3',
                        'te': 'keV',
                        'ti': 'keV',
                        'vr': 'cm/s',
                        'vt': 'cm/s',
                        'vz': 'cm/s',
                        'r': 'cm',
                        'z': 'cm',
                        'r2d': 'cm',
                        'z2d': 'cm'}

        # Define dimension names for xarray/pandas compatibility
        plasma_dim_names = {'dene': ['r', 'z'],
                           'denimp': ['r', 'z'],
                           'deni': ['species', 'r', 'z'],
                           'te': ['r', 'z'],
                           'ti': ['r', 'z'],
                           'zeff': ['r', 'z'],
                           'denn': ['r', 'z'],
                           'vr': ['r', 'z'],
                           'vt': ['r', 'z'],
                           'vz': ['r', 'z'],
                           'r2d': ['r', 'z'],
                           'z2d': ['r', 'z'],
                           'mask': ['r', 'z'],
                           'species_mass': ['species']}

        # Mark coordinate arrays as dimension scales
        # Include species_mass as the scale for the species dimension
        plasma_coord_scales = {'r': 'r', 'z': 'z', 'species_mass': 'species'}

        # Add dimension names for profiles subgroup if it might exist
        # These will use the nested key format that write_data now supports
        plasma_dim_names.update({
            'profiles/dene': ['rho'],
            'profiles/te': ['rho'],
            'profiles/ti': ['rho'],
            'profiles/zeff': ['rho'],
            'profiles/omega': ['rho'],
            'profiles/rho': ['rho']
        })

        # Mark rho as a coordinate scale for the profiles
        plasma_coord_scales['profiles/rho'] = 'rho'

        # Add descriptions for profiles subgroup (if it exists)
        plasma_description.update({
            'profiles/dene': 'Electron density profile',
            'profiles/te': 'Electron temperature profile',
            'profiles/ti': 'Ion temperature profile',
            'profiles/zeff': 'Effective charge profile',
            'profiles/omega': 'Toroidal rotation profile',
            'profiles/rho': 'Normalized radial coordinate'
        })

        # Add units for profiles subgroup (if it exists)
        plasma_units.update({
            'profiles/dene': 'cm^-3',
            'profiles/te': 'keV',
            'profiles/ti': 'keV',
            'profiles/omega': 'rad/s'
        })

        write_data(g_plasma, plasma, desc = plasma_description, units = plasma_units,
                   dim_names=plasma_dim_names, coord_scales=plasma_coord_scales, name='plasma')

        # Create fields group
        g_fields = hf.create_group('fields')

        # Electromagnetic fields attributes
        g_fields.attrs['description'] = 'Electromagnetic Fields'
        g_fields.attrs['coordinate_system'] = 'Cylindrical'

        fields_description = {'data_source': 'Source of the EM equilibrium',
                              'mask': 'Boolean mask that indicates where the fields are well defined',
                              'time': 'Time',
                              'br': 'Magnetic field in the r-direction: Br(r,z)',
                              'bt': 'Magnetic field in the theta/torodial-direction: Bt(r,z)',
                              'bz': 'Magnetic field in the z-direction: Bz(r,z)',
                              'er': 'Electric field in the r-direction: Er(r,z)',
                              'et': 'Electric field in the theta/torodial-direction: Et(r,z)',
                              'ez': 'Electric field in the z-direction: Ez(r,z)',
                              'nr': 'Number of R values',
                              'nz': 'Number of Z values',
                              'r': 'Radius',
                              'z': 'Z',
                              'r2d': 'Radius grid: R(r,z)',
                              'z2d': 'Z grid: Z(r,z)'}

        fields_units = {'time': 's',
                        'br': 'T',
                        'bt': 'T',
                        'bz': 'T',
                        'er': 'V/m',
                        'et': 'V/m',
                        'ez': 'V/m',
                        'nr': 'V/m',
                        'nz': 'V/m',
                        'r': 'cm',
                        'z': 'cm',
                        'r2d': 'cm',
                        'z2d': 'cm'}

        # Define dimension names for xarray/pandas compatibility
        fields_dim_names = {'br': ['r', 'z'],
                           'bt': ['r', 'z'],
                           'bz': ['r', 'z'],
                           'er': ['r', 'z'],
                           'et': ['r', 'z'],
                           'ez': ['r', 'z'],
                           'r2d': ['r', 'z'],
                           'z2d': ['r', 'z'],
                           'mask': ['r', 'z']}  # Added mask to ensure consistency

        # Mark coordinate arrays as dimension scales
        fields_coord_scales = {'r': 'r', 'z': 'z'}

        write_data(g_fields, fields, desc = fields_description, units = fields_units,
                   dim_names=fields_dim_names, coord_scales=fields_coord_scales, name='fields')

    if os.path.isfile(filename):
        success('Equilibrium file created: '+filename)
    else:
        error('Equilibrium file creation failed.')

def write_distribution(filename, distri):
    """
    #+#write_distribution
    #+Write fast-ion distribution to a HDF5 file
    #+***
    #+##Input Arguments
    #+     **filename**: Name of the distribution file
    #+
    #+     **dist**: Fast-ion distribution distionary
    #+
    #+##Example Usage
    #+```dist
    #+>>> write_distribution(filename, distri)
    #+```
    """
    info('Writing fast-ion distribution file...')

    description = {'data_source': 'Source of the fast-ion distribution',
                   'type': 'Distribution type: 1="Guiding Center Density Function", 2="Guiding Center ' \
                           'Monte Carlo", 3="Full Orbit Monte Carlo"',
                   'time': 'Distribution time'}

    units = {'time': 's'}

    if distri['type'] == 1:
        description['nenergy'] = 'Number of energy values'
        description['npitch'] = 'Number of pitch values'
        description['energy'] = 'Energy'
        description['pitch'] = 'Pitch: p = v_parallel/v  w.r.t. the magnetic field'
        description['f'] = 'Fast-ion density function: F(E,p,R,Z)'
        description['denf'] = 'Fast-ion density: Denf(r,z)'
        description['nr'] = 'Number of R values'
        description['nz'] = 'Number of Z values'
        description['r'] = 'Radius'
        description['z'] = 'Z'
        description['r2d'] = 'Radius grid: R(r,z)'
        description['z2d'] = 'Z grid: Z(r,z)'

        units['energy'] = 'keV'
        units['f'] = 'fast-ions/(dE*dP*cm^3)'
        units['denf'] = 'cm^-3'
        units['r'] = 'cm'
        units['z'] = 'cm'
        units['r2d'] = 'cm'
        units['z2d'] = 'cm'
    else:
        description['nparticle'] = 'Number of MC particles'
        description['nclass'] = 'Number of orbit classes'
        description['r'] = 'R position of a MC particle'
        description['z'] = 'Z position of a MC particle'
        description['weight'] = 'Weight of a MC particle: sum(weight) = # of fast-ions '
        description['class'] = 'Orbit class of a MC particle: class in Set(1:nclass)'

        units['r'] = 'cm'
        units['z'] = 'cm'
        units['weight'] = 'fast-ions/particle'

        if distri['type'] == 2:
            description['energy'] = 'Energy of a MC particle'
            description['pitch'] = 'Pitch of a MC particle: p = v_parallel/v  w.r.t. the magnetic field'
        else:
            description['vr'] ='Radial velocity of a MC particle'
            description['vt'] = 'Torodial velocity of a MC particle'
            description['vz'] = 'Z velocity of a MC particle'

            units['vr'] = 'cm/s'
            units['vt'] = 'cm/s'
            units['vz'] = 'cm/s'

    with h5py.File(filename, 'w') as hf:
        # File attr
        hf.attrs['description'] = 'Fast-ion distribution for FIDASIM'
        hf.attrs['coordinate_system'] = 'Cylindrical'

        # Define dimension names for xarray/pandas compatibility
        dim_names = {}
        coord_scales = {}

        if distri['type'] == 1:
            # Guiding Center Density Function
            dim_names['f'] = ['energy', 'pitch', 'r', 'z']
            dim_names['denf'] = ['r', 'z']
            dim_names['r2d'] = ['r', 'z']
            dim_names['z2d'] = ['r', 'z']
            coord_scales = {'r': 'r', 'z': 'z', 'energy': 'energy', 'pitch': 'pitch'}
        else:
            # Monte Carlo distributions
            dim_names['r'] = ['particle']
            dim_names['z'] = ['particle']
            dim_names['weight'] = ['particle']
            dim_names['class'] = ['particle']
            if distri['type'] == 2:
                dim_names['energy'] = ['particle']
                dim_names['pitch'] = ['particle']
            else:
                dim_names['vr'] = ['particle']
                dim_names['vt'] = ['particle']
                dim_names['vz'] = ['particle']

        write_data(hf, distri, desc = description, units=units,
                   dim_names=dim_names, coord_scales=coord_scales, name='distribution')

    if os.path.isfile(filename):
        success('Distribution file created: ' + filename)
    else:
        error('Distribution file creation failed.')

def prefida(inputs, grid, nbi, plasma, fields, fbm, spec=None, npa=None, nc=None, cfpd=None, use_abs_path=True):
    """
    #+#prefida
    #+Checks FIDASIM inputs and writes FIDASIM input files
    #+***
    #+##Input Arguments
    #+     **inputs**: Inputs structure
    #+
    #+     **grid**: Interpolation grid structure
    #+
    #+     **nbi**: Neutral beam geometry structure (can be None or empty dict for passive-only runs)
    #+
    #+     **plasma**: Plasma parameters structure
    #+
    #+     **fields**: Electromagnetic fields structure
    #+
    #+     **fbm**: Fast-ion distribution structure
    #+
    #+##Keyword Arguments
    #+     **spec**: Optional, Spectral geometry structure
    #+
    #+     **npa**: Optional, NPA geometry structure
    #+
    #+     **nc**: Optional, NC geometry structure
    #+
    #+     **cfpd**: Optional, CFPD geometry structure
    #+
    #+##Example Usage
    #+```python
    #+>>> # With beam (active measurements)
    #+>>> prefida(inputs, grid, nbi, plasma, fields, fbm, spec=spec, npa=npa, nc=nc)
    #+>>> # Passive-only (no beam) - pass None or empty dict
    #+>>> prefida(inputs, grid, None, plasma, fields, fbm, spec=spec, npa=npa, nc=nc)
    #+```
    """
    # CHECK INPUTS
    inputs = check_inputs(inputs, use_abs_path=use_abs_path)

    # MAKE DIRECTORIES IF THEY DONT EXIST
    if not os.path.isdir(inputs['result_dir']):
        os.makedirs(inputs['result_dir'])

    # CHECK INTERPOLATION GRID
    check_grid(grid)

    # CHECK BEAM INPUTS (optional for passive-only runs)
    # Accept None or empty dict for passive-only mode
    if nbi is not None and nbi:  # Check if nbi exists and is not empty
        nbi = check_beam(inputs, nbi)
    else:
        # Passive-only mode: nbi is None or empty dict
        # Check if beam parameters are present in inputs
        has_beam_params = all(key in inputs for key in ['nx', 'ny', 'nz', 'einj', 'pinj', 'current_fractions'])
        if not has_beam_params:
            warn('='*70)
            warn('PASSIVE-ONLY MODE DETECTED')
            warn('No beam geometry provided (nbi=None or {}) and no beam parameters in inputs.')
            warn('Only passive diagnostics will be calculated:')
            warn('  - Passive FIDA (calc_pfida)')
            warn('  - Passive NPA (calc_pnpa)')
            warn('  - Neutron rates (calc_neutron)')
            warn('Active diagnostics (FIDA, NPA, BES, DCX, HALO) will NOT be available.')
            warn('='*70)

    # CHECK PLASMA PARAMETERS
    plasma = check_plasma(inputs, grid, plasma)

    # CHECK ELECTROMAGNETIC FIELDS
    fields = check_fields(inputs, grid, fields)

    # CHECK FAST-ION DISTRIBUTION
    fbm = check_distribution(inputs, grid, fbm)

    # CHECK FIDA/BES
    if spec is not None:
        check_spec(inputs, spec)

    # CHECK NPA
    if npa is not None:
        check_npa(inputs, npa)
    
    # CHECK NC
    if nc is not None:
        check_nc(inputs, nc)

    # CHECK CFPD
    if cfpd is not None:
        check_cfpd(inputs, cfpd)

    # WRITE FIDASIM INPUT FILES
    write_namelist(inputs['input_file'], inputs)

    # WRITE GEOMETRY FILE
    write_geometry(inputs['geometry_file'], nbi, spec=spec, npa=npa, nc=nc, cfpd=cfpd)

    # WRITE EQUILIBRIUM FILE
    write_equilibrium(inputs['equilibrium_file'], plasma, fields)

    # WRITE DISTRIBUTION FILE
    write_distribution(inputs['distribution_file'], fbm)

    print('')
    print('')
    success('FIDASIM pre-processing completed')
    print('To run FIDASIM use the following command')
    print(get_fidasim_dir() + os.sep + 'fidasim ' + inputs['result_dir'] + os.sep + inputs['runid'] + '_inputs.dat')
    print('')
    print('')
