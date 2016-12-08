#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import os
from lib.python_prefida.success import success
from lib.python_prefida.info import info
from lib.python_prefida.error import error
from lib.python_prefida.write_data import write_data


def write_geometry(filename, nbi, spec=None, npa=None):
    #+#write_geometry
    #+Write geometry values to a HDF5 file
    #+***
    #+##Input Arguments
    #+     **filename**: Name of the geometry file
    #+
    #+     **nbi**: NBI geometry structure
    #+
    #+##Keyword Arguments
    #+     **spec**: Optional, Spectral geometry structure
    #+
    #+     **npa**: Optional, NPA geometry structure
    #+
    #+##Example Usage
    #+```idl
    #+IDL> write_geometry, filename, nbi, spec=spec, npa=npa
    #+```
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

        write_data(g_nbi, nbi, nbi_description, nbi_units, name='nbi')

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

            write_data(g_spec, spec, spec_description, spec_units, name='spec')

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

            write_data(g_npa, npa, npa_description, npa_units, name='npa')

    if os.path.isfile(filename):
        success('Geometry file created: ' + filename)
    else:
        error('Geometry file creation failed.')
