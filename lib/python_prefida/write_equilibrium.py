#!/usr/bin/env python
# -*- coding: utf-8 -*-

from info import info
import h5py


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
                              'time': 'Time'
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


    plasma_units =
    plasma_time_unit = {attribute,obj:'/plasma/time', $
                        name:'units', $
                        data:'s'}

    plasma_dene_unit = {attribute,obj:'/plasma/dene', $
                        name:'units', $
                        data:'cm^-3'}
    plasma_te_unit = {attribute,obj:'/plasma/te', $
                      name:'units', $
                      data:'keV'}
    plasma_ti_unit = {attribute,obj:'/plasma/ti', $
                      name:'units', $
                      data:'keV'}
    plasma_vr_unit = {attribute,obj:'/plasma/vr', $
                      name:'units', $
                      data:'cm/s'}
    plasma_vt_unit = {attribute,obj:'/plasma/vt', $
                      name:'units', $
                      data:'cm/s'}
    plasma_vz_unit = {attribute,obj:'/plasma/vz', $
                      name:'units', $
                      data:'cm/s'}
    plasma_r_unit = {attribute,obj:'/plasma/r', $
                     name:'units', $
                     data:'cm'}
                      plasma_z_unit = {attribute,obj:'/plasma/z', $
                     name:'units', $
                     data:'cm'}
                     plasma_r2d_unit = {attribute,obj:'/plasma/r2d', $
                       name:'units', $
                       data:'cm'}
                      plasma_z2d_unit = {attribute,obj:'/plasma/z2d', $
                       name:'units', $
                       data:'cm'}



    ;; Electromagnetic fields attributes
    fields_desc = {attribute,obj:'/fields', $
                   name:'description', $
                   data:'Electromagnetic Fields'}

    fields_cs = {attribute,obj:'/fields', $
                 name:'coordinate_system', $
                 data:'Cylindrical'}

    fields_ds_desc = {attribute,obj:'/fields/data_source', $
                      name:'description', $
                      data:'Source of the EM equilibrium'}

    fields_mask_desc = {attribute,obj:'/fields/mask',$
                        name:'description', $
                        data:'Boolean mask that indicates where' +$
                        ' the fields are well defined'}

    fields_time_desc = {attribute,obj:'/fields/time', $
                        name:'description', $
                        data:'Time'}
    fields_time_unit = {attribute,obj:'/fields/time', $
                        name:'units', $
                        data:'s'}

    fields_br_desc = {attribute,obj:'/fields/br', $
                      name:'description', $
                      data:'Magnetic field in the r-direction: Br(r,z)'}
    fields_br_unit = {attribute,obj:'/fields/br', $
                      name:'units', $
                      data:'T'}

    fields_bt_desc = {attribute,obj:'/fields/bt', $
                      name:'description', $
                      data:'Magnetic field in the theta/torodial-direction: Bt(r,z)'}
    fields_bt_unit = {attribute,obj:'/fields/bt', $
                      name:'units', $
                      data:'T'}

    fields_bz_desc = {attribute,obj:'/fields/bz', $
                      name:'description', $
                      data:'Magnetic field in the z-direction: Bz(r,z)'}
    fields_bz_unit = {attribute,obj:'/fields/bz', $
                      name:'units', $
                      data:'T'}

    fields_er_desc = {attribute,obj:'/fields/er', $
                      name:'description', $
                      data:'Electric field in the r-direction: Er(r,z)'}
    fields_er_unit = {attribute,obj:'/fields/er', $
                      name:'units', $
                      data:'V/m'}

    fields_et_desc = {attribute,obj:'/fields/et', $
                      name:'description', $
                      data:'Electric field in the theta/torodial-direction: Et(r,z)'}
    fields_et_unit = {attribute,obj:'/fields/et', $
                      name:'units', $
                      data:'V/m'}

    fields_ez_desc = {attribute,obj:'/fields/ez', $
                      name:'description', $
                      data:'Electric field in the z-direction: Ez(r,z)'}
    fields_ez_unit = {attribute,obj:'/fields/ez', $
                      name:'units', $
                      data:'V/m'}

    fields_nr_desc = {attribute,obj:'/fields/nr', $
                      name:'description', $
                      data:'Number of R values'}

    fields_nz_desc = {attribute,obj:'/fields/nz', $
                      name:'description', $
                      data:'Number of Z values'}

    fields_r_desc = {attribute,obj:'/fields/r', $
                     name:'description', $
                     data:'Radius'}
    fields_r_unit = {attribute,obj:'/fields/r', $
                     name:'units', $
                     data:'cm'}

    fields_z_desc = {attribute,obj:'/fields/z', $
                     name:'description', $
                     data:'Z'}
    fields_z_unit = {attribute,obj:'/fields/z', $
                     name:'units', $
                     data:'cm'}

    fields_r2d_desc = {attribute,obj:'/fields/r2d', $
                       name:'description', $
                       data:'Radius grid: R(r,z)'}
    fields_r2d_unit = {attribute,obj:'/fields/r2d', $
                       name:'units', $
                       data:'cm'}

    fields_z2d_desc = {attribute,obj:'/fields/z2d', $
                       name:'description', $
                       data:'Z grid: Z(r,z)'}
    fields_z2d_unit = {attribute,obj:'/fields/z2d', $
                       name:'units', $
                       data:'cm'}

    fields_atts = [fields_desc, fields_cs, fields_ds_desc, $
                   fields_mask_desc, $
                   fields_time_desc, fields_time_unit, $
                   fields_br_desc, fields_br_unit, $
                   fields_bt_desc, fields_bt_unit, $
                   fields_bz_desc, fields_bz_unit, $
                   fields_er_desc, fields_er_unit, $
                   fields_et_desc, fields_et_unit, $
                   fields_ez_desc, fields_ez_unit, $
                   fields_nr_desc, fields_nz_desc, $
                   fields_r_desc, fields_r_unit, $
                   fields_z_desc, fields_z_unit, $
                   fields_r2d_desc, fields_r2d_unit, $
                   fields_z2d_desc, fields_z2d_unit ]

    atts = [root_atts, plasma_atts, fields_atts]

    write_hdf5,["plasma","fields"],filename=filename, atts=atts, /clobber

    if file_test(filename) then begin
        success, 'Equilibrium file created: '+filename
    endif else begin
        error, 'Equilibrium file creation failed.'
    endelse
END
