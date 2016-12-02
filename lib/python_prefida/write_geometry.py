#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import os
from success import success
from info import info
from error import error


def write_geometry(filename, nbi, spec=None, npa=None):
    """Brief Description

    Expanded description

    Sample usage:
    -------------
    >>> Description

    Notes
    -----

    Parameters
    ----------
    parameter : type

        Description

    keyword : type

        Description

    Returns
    -------
    result : type

        Description

    History
    -------
    Created on Fri Dec  2 09:44:20 2016 by nbolte

    To Do
    -----
    * Consider setting h5 driver according to current OS

    """
    """
    #write_geometry
    Write geometry values to a HDF5 file with up to 3 groups: nbi, spec, npa
    ***
    ##Input Arguments
         **filename**: Name of the geometry file

         **nbi**: NBI geometry structure

    ##Keyword Arguments
         **spec**: Optional, Spectral geometry structure

         **npa**: Optional, NPA geometry structure

    ##Example Usage
    ```idl
    IDL> write_geometry, filename, nbi, spec=spec, npa=npa
    ```
    """
    info('Writing geometry file...')

    # Create and open h5 file
    with h5py.File(filename, 'w') as hf:
        # File attributes
        hf.attrs['description'] = 'Geometric quantities for FIDASIM'

        # Create nbi group
        g_nbi = hf.create_group('nbi')

        # nbi att
        g_nbi.attrs['description'] = 'Neutral Beam Geometry'
        g_nbi.attrs['coordinate_system'] = 'Right-handed cartesian'

        # Save nbi data
        key = 'data_source'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Source of the NBI geometry'

        key = 'name'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Beam name'

        key = 'src'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Position of the center of the beam source grid'
        ds.attrs['units'] = 'cm'

        key = 'axis'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Axis of the beam centerline: Centerline(t) = src + axis*t '
        ds.attrs['units'] = 'cm'

        key = 'focy'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Horizonal focal length of the beam'
        ds.attrs['units'] = 'cm'

        key = 'focz'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Vertical focal length of the beam'
        ds.attrs['units'] = 'cm'

        key = 'divy'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Horizonal divergences of the beam. One for each energy component'
        ds.attrs['units'] = 'radians'

        key = 'divz'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Vertical divergences of the beam. One for each energy component'
        ds.attrs['units'] = 'radians'

        key = 'widy'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Half width of the beam source grid'
        ds.attrs['units'] = 'cm'

        key = 'widz'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Half height of the beam source grid'
        ds.attrs['units'] = 'cm'

        key = 'shape'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Shape of the beam source grid: 1="rectangular", 2="circular"'

        key = 'naperture'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Number of apertures'

        key = 'ashape'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Shape of the aperture(s): 1="rectangular", 2="circular"'

        key = 'awidy'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Half width of the aperture(s)'
        ds.attrs['units'] = 'cm'

        key = 'awidz'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Half height of the aperture(s)'
        ds.attrs['units'] = 'cm'

        key = 'aoffy'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Horizontal (y) offset of the aperture(s) relative to the +x aligned beam centerline'
        ds.attrs['units'] = 'cm'

        key = 'aoffz'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Vertical (z) offset of the aperture(s) relative to the +x aligned beam centerline'
        ds.attrs['units'] = 'cm'

        key = 'adist'
        ds = g_nbi.create_dataset(key, data = nbi[key])
        ds.attrs['description'] = 'Distance from the center of the beam source grid to the aperture(s) plane'
        ds.attrs['units'] = 'cm'

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

            for key in spec:
                # Create dataset
                ds = g_spec.create_dataset(key, data = spec[key])

                # Add descrption attr
                ds.attrs['description'] = spec_description[key]

                # Add units attr
                if key in spec_units:
                    ds.attrs['units'] = spec_units[key]

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

            for key in npa:
                # Create dataset
                ds = g_npa.create_dataset(key, data = npa[key])

                # Add descrption attr
                ds.attrs['description'] = npa_description[key]

                # Add units attr
                if key in npa_units:
                    ds.attrs['units'] = npa_units[key]

    if os.path.isfile(filename):
        success('Geometry file created: ' + filename)
    else:
        error('Geometry file creation failed.')
