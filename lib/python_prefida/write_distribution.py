#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import os
from lib.python_prefida.success import success
from lib.python_prefida.error import error
from lib.python_prefida.info import info


def write_distribution(filename, distri):
    #+##`write_distribution, filename, dist`
    #+Write fast-ion distribution to a HDF5 file
    #+
    #+###Input Arguments
    #+     **filename**: Name of the distribution file
    #+
    #+     **dist**: Fast-ion distribution structure
    #+
    #+###Example Usage
    #+```idl
    #+IDL> write_distribution, filename, distri
    #+```
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

        for key in distri:
            # Create dataset
            ds = hf.create_dataset(key, data = distri[key])

            # Add descrption attr
            ds.attrs['description'] = description[key]

            # Add units attr
            if key in units:
                ds.attrs['units'] = units[key]

    if os.path.isfile(filename):
        success('Distribution file created: ' + filename)
    else:
        error('Distribution file creation failed.')
