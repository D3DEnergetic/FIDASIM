#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib.python_prefida import info
import h5py
import os
from lib.python_prefida import success
from lib.python_prefida import error


def write_equilibrium(filename, plasma, fields):
    """Brief Description

    Expanded description

    Sample usage:
    -------------
    >>> Description

    Notes
    -----
    Two groups: 'plasma' and 'fields'

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
    Created on Fri Dec  2 13:57:20 2016 by nbolte

    To Do
    -----

    ;+##`write_equilibrium, filename, plasma, fields`
    ;+Write MHD equilibrium values to a HDF5 file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the equilibrium file
    ;+
    ;+     **plasma**: Plasma structure
    ;+
    ;+     **fields**: Electromagnetic fields structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_equilibrium, filename, plasma, fields
    ;+```
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
                              'dene': 'Electron Number Density: Dene(r,z)',
                              'te': 'Electron Temperature: Te(r,z)',
                              'ti': 'Ion Temperature: Ti(r,z)',
                              'zeff': 'Effective Nuclear Charge: Zeff(r,z)',
                              'vr': 'Bulk plasma flow in the r-direction: Vr(r,z)',
                              'vt': 'Bulk plasma flow in the theta/torodial-direction: Vt(r,z)',
                              'vz': 'Bulk plasma flow in the z-direction: Vz(r,z)',
                              'nr': 'Number of R values',
                              'nz': 'Number of Z values',
                              'r': 'Radius',
                              'z': 'Z',
                              'r2d': 'Radius grid: R(r,z)',
                              'z2d': 'Z grid: Z(r,z)',
                              'mask': 'Boolean mask that indicates where the plasma parameters are well defined'}

        plasma_units = {'time': 's',
                        'dene': 'cm^-3',
                        'te': 'keV',
                        'ti': 'keV',
                        'vr': 'cm/s',
                        'vt': 'cm/s',
                        'vz': 'cm/s',
                        'r': 'cm',
                        'z': 'cm',
                        'r2d': 'cm',
                        'z2d': 'cm'}

        for key in plasma:
                # Create dataset
                ds = g_plasma.create_dataset(key, data = plasma[key])

                # Add descrption attr
                ds.attrs['description'] = plasma_description[key]

                # Add units attr
                if key in plasma_units:
                    ds.attrs['units'] = plasma_units[key]

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

    for key in fields:
                # Create dataset
                ds = g_fields.create_dataset(key, data = fields[key])

                # Add descrption attr
                ds.attrs['description'] = fields_description[key]

                # Add units attr
                if key in fields_units:
                    ds.attrs['units'] = fields_units[key]

    if os.path.isfile(filename):
        success('Equilibrium file created: '+filename)
    else:
        error('Equilibrium file creation failed.')