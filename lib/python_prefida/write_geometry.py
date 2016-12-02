#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py


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
        # Root att
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

        nbi_shape_desc = {attribute,obj:'/nbi/shape',
                          name:'description',
                          data:'Shape of the beam source grid: 1="rectangular", 2="circular"'}

        nbi_naperture_desc = {attribute,obj:'/nbi/naperture',
                              name:'description',
                              data:'Number of apertures'}

        nbi_ashape_desc = {attribute,obj:'/nbi/ashape',
                           name:'description',
                           data:'Shape of the aperture(s): 1="rectangular", 2="circular"'}

        nbi_awidy_desc = {attribute,obj:'/nbi/awidy',
                          name:'description',
                          data:'Half width of the aperture(s)'}
        nbi_awidy_unit = {attribute,obj:'/nbi/awidy',
                          name:'units',
                          data:'cm'}

        nbi_awidz_desc = {attribute,obj:'/nbi/awidz',
                          name:'description',
                          data:'Half height of the aperture(s)'}
        nbi_awidz_unit = {attribute,obj:'/nbi/awidz',
                          name:'units',
                          data:'cm'}

        nbi_aoffy_desc = {attribute,obj:'/nbi/aoffy',
                          name:'description',
                          data:'Horizontal (y) offset of the aperture(s) relative to the +x aligned beam centerline'}
        nbi_aoffy_unit = {attribute,obj:'/nbi/aoffy',
                          name:'units',
                          data:'cm'}

        nbi_aoffz_desc = {attribute,obj:'/nbi/aoffz',
                          name:'description',
                          data:'Vertical (z) offset of the aperture(s) relative to the +x aligned beam centerline'}
        nbi_aoffz_unit = {attribute,obj:'/nbi/aoffz',
                          name:'units',
                          data:'cm'}

        nbi_adist_desc = {attribute,obj:'/nbi/adist',
                          name:'description',
                          data:'Distance from the center of the beam source grid to the aperture(s) plane'}
        nbi_adist_unit = {attribute,obj:'/nbi/adist',
                          name:'units',
                          data:'cm'}

    nbi_atts = [nbi_desc, nbi_cs, nbi_ds_desc,nbi_name_desc,
                nbi_shape_desc, nbi_src_desc, nbi_src_unit,
                nbi_axis_desc, nbi_axis_unit,
                nbi_focy_desc, nbi_focy_unit,
                nbi_focz_desc, nbi_focz_unit,
                nbi_divy_desc, nbi_divy_unit,
                nbi_divz_desc, nbi_divz_unit,
                nbi_widy_desc, nbi_widy_unit,
                nbi_widz_desc, nbi_widz_unit,
                nbi_naperture_desc, nbi_ashape_desc,
                nbi_awidy_desc, nbi_awidy_unit,
                nbi_awidz_desc, nbi_awidz_unit,
                nbi_aoffy_desc, nbi_aoffy_unit,
                nbi_aoffz_desc, nbi_aoffz_unit,
                nbi_adist_desc, nbi_adist_unit ]

    if spec is not None:
            # Create spec group
            g_spec = hf.create_group('spec')

            # Save spec data
            for key in spec:
                g_spec.create_dataset('spec/' + key, data = spec[key])

        if npa is not None:
            # Create npa group
            g_npa = hf.create_group('npa')

            # Save npa data
            for key in npa:
                g_npa.create_dataset('npa/' + key, data = npa[key])

    # Spectroscopic attributes
    spec_desc = {attribute,obj:'/spec',
                 name:'description',
                 data:'FIDA/BES Chord Geometry'}
    spec_cs = {attribute,obj:'/spec',
               name:'coordinate_system',
               data:'Right-handed cartesian'}

    spec_ds_desc = {attribute,obj:'/spec/data_source',
                    name:'description',
                    data:'Source of the chord geometry'}

    spec_nchan_desc = {attribute,obj:'/spec/nchan',
                       name:'description',
                       data:'Number of channels'}

    spec_system_desc = {attribute,obj:'/spec/system',
                        name:'description',
                        data:'Names of the different spectrocopic systems'}

    spec_id_desc = {attribute,obj:'/spec/id',
                    name:'description',
                    data:'Line of sight ID'}

    spec_lens_desc = {attribute,obj:'/spec/lens',
                      name:'description',
                      data:'Positions of the lenses'}
    spec_lens_unit = {attribute,obj:'/spec/lens',
                      name:'units',
                      data:'cm'}

    spec_axis_desc = {attribute,obj:'/spec/axis',
                      name:'description',
                      data:'Optical axis of the lines of sight: LOS(t) = lens + axis*t '}
    spec_axis_unit = {attribute,obj:'/spec/axis',
                      name:'units',
                      data:'cm'}

    spec_radius_desc = {attribute,obj:'/spec/radius',
                        name:'description',
                        data:'Line of sight radius at midplane or tangency point' }
    spec_radius_unit = {attribute,obj:'/spec/radius',
                        name:'units',
                        data:'cm'}

    spec_sigma_desc = {attribute,obj:'/spec/sigma_pi',
                       name:'description',
                       data:'Ratio of the intensities of the sigma and pi stark lines. Measured quantity'}

    spec_spot_desc = {attribute,obj:'/spec/spot_size',
                      name:'description',
                      data:'Radius of spot size'}
    spec_spot_unit = {attribute,obj:'/spec/spot_size',
                      name:'units',
                      data:'cm'}

    spec_atts = [spec_desc,spec_cs, spec_ds_desc,
                 spec_nchan_desc, spec_system_desc,
                 spec_id_desc,
                 spec_lens_desc,spec_lens_unit,
                 spec_axis_desc,spec_axis_unit,
                 spec_radius_desc, spec_radius_unit,
                 spec_sigma_desc,
                 spec_spot_desc, spec_spot_unit]

    #NPA attributes
    npa_desc = {attribute,obj:'/npa',
                name:'description',
                data:'NPA Geometry'}

    npa_cs = {attribute,obj:'/npa',
              name:'coordinate_system',
              data:'Right-handed cartesian'}

    npa_ds_desc = {attribute,obj:'/npa/data_source',
                   name:'description',
                   data:'Source of the NPA geometry'}

    npa_nchan_desc = {attribute,obj:'/npa/nchan',
                      name:'description',
                      data:'Number of channels'}

    npa_system_desc = {attribute,obj:'/npa/system',
                       name:'description',
                       data:'Names of the different NPA systems'}

    npa_id_desc = {attribute,obj:'/npa/id',
                   name:'description',
                   data:'Line of sight ID'}

    npa_dshape_desc = {attribute,obj:'/npa/d_shape',
                       name:'description',
                       data:'Shape of the detector: 1="rectangular", 2="circular"'}

    npa_dcent_desc = {attribute,obj:'/npa/d_cent',
                      name:'description',
                      data:'Center of the detector'}
    npa_dcent_unit = {attribute,obj:'/npa/d_cent',
                      name:'units',
                      data:'cm'}

    npa_dtedge_desc = {attribute,obj:'/npa/d_tedge',
                       name:'description',
                       data:'Center of the detectors top edge'}
    npa_dtedge_unit = {attribute,obj:'/npa/d_tedge',
                       name:'units',
                       data:'cm'}

    npa_dredge_desc = {attribute,obj:'/npa/d_redge',
                       name:'description',
                       data:'Center of the detectors right edge'}
    npa_dredge_unit = {attribute,obj:'/npa/d_redge',
                       name:'units',
                       data:'cm'}

    npa_ashape_desc = {attribute,obj:'/npa/a_shape',
                       name:'description',
                       data:'Shape of the aperture: 1="rectangular", 2="circular"'}

    npa_acent_desc = {attribute,obj:'/npa/a_cent',
                      name:'description',
                      data:'Center of the aperture'}
    npa_acent_unit = {attribute,obj:'/npa/a_cent',
                      name:'units',
                      data:'cm'}

    npa_atedge_desc = {attribute,obj:'/npa/a_tedge',
                       name:'description',
                       data:'Center of the apertures top edge'}
    npa_atedge_unit = {attribute,obj:'/npa/a_tedge',
                       name:'units',
                       data:'cm'}

    npa_aredge_desc = {attribute,obj:'/npa/a_redge',
                       name:'description',
                       data:'Center of the apertures right edge'}
    npa_aredge_unit = {attribute,obj:'/npa/a_redge',
                       name:'units',
                       data:'cm'}

    npa_radius_desc = {attribute,obj:'/npa/radius',
                       name:'description',
                       data:'Line of sight radius at midplane or tangency point' }
    npa_radius_unit = {attribute,obj:'/npa/radius',
                       name:'units',
                       data:'cm'}

    npa_atts = [npa_desc, npa_cs, npa_ds_desc,
                npa_nchan_desc, npa_system_desc,
                npa_id_desc,
                npa_dshape_desc, npa_ashape_desc,
                npa_dcent_desc, npa_dcent_unit,
                npa_acent_desc, npa_acent_unit,
                npa_dtedge_desc, npa_dtedge_unit,
                npa_atedge_desc, npa_atedge_unit,
                npa_dredge_desc, npa_dredge_unit,
                npa_aredge_desc, npa_aredge_unit,
                npa_radius_desc, npa_radius_unit ]

    atts = [root_atts, nbi_atts]
    geom = {nbi:nbi}

    if spec is not None:
        geom['spec'] = spec
        atts = [atts, spec_atts]


    if npa is not None:
        geom['npa'] = npa
        atts = [atts,npa_atts]


    write_hdf5, geom, filename=filename, atts=atts, /clobber

    if file_test(filename):
        success('Geometry file created: ' + filename)
    else:
        error('Geometry file creation failed.')




def write_equilibrium(filename, plasma, fields):
    """
    ##`write_equilibrium, filename, plasma, fields`
    Write MHD equilibrium values to a HDF5 file

    ###Input Arguments
         **filename**: Name of the equilibrium file

         **plasma**: Plasma structure

         **fields**: Electromagnetic fields structure

    ###Example Usage
    ```idl
    IDL> write_equilibrium, filename, plasma, fields
    ```
    """
    info('Writing equilibrium file...')

    root_atts = {attribute,obj:'/',
                 name:'description',
                 data:'Plasma Parameters and Electromagnetic Fields for FIDASIM'}

    # Plasma Attributes
    plasma_desc = {attribute,obj:'/plasma',
                   name:'description',
                   data:'Plasma Parameters'}

    plasma_cs = {attribute,obj:'/plasma',
                 name:'coordinate_system',
                 data:'Cylindrical'}

    plasma_ds_desc = {attribute, obj:'/plasma/data_source',
                      name:'description',
                      data:'Source of the plasma parameters'}

    plasma_time_desc = {attribute,obj:'/plasma/time',
                        name:'description',
                        data:'Time'}
    plasma_time_unit = {attribute,obj:'/plasma/time',
                        name:'units',
                        data:'s'}

    plasma_dene_desc = {attribute,obj:'/plasma/dene',
                        name:'description',
                        data:'Electron Number Density: Dene(r,z)'}
    plasma_dene_unit = {attribute,obj:'/plasma/dene',
                        name:'units',
                        data:'cm^-3'}

    plasma_te_desc = {attribute,obj:'/plasma/te',
                      name:'description',
                      data:'Electron Temperature: Te(r,z)'}
    plasma_te_unit = {attribute,obj:'/plasma/te',
                      name:'units',
                      data:'keV'}

    plasma_ti_desc = {attribute,obj:'/plasma/ti',
                      name:'description',
                      data:'Ion Temperature: Ti(r,z)'}
    plasma_ti_unit = {attribute,obj:'/plasma/ti',
                      name:'units',
                      data:'keV'}

    plasma_zeff_desc = {attribute,obj:'/plasma/zeff',
                        name:'description',
                        data:'Effective Nuclear Charge: Zeff(r,z)'}

    plasma_vr_desc = {attribute,obj:'/plasma/vr',
                      name:'description',
                      data:'Bulk plasma flow in the r-direction: Vr(r,z)'}
    plasma_vr_unit = {attribute,obj:'/plasma/vr',
                      name:'units',
                      data:'cm/s'}

    plasma_vt_desc = {attribute,obj:'/plasma/vt',
                      name:'description',
                      data:'Bulk plasma flow in the theta/torodial-direction: Vt(r,z)'}
    plasma_vt_unit = {attribute,obj:'/plasma/vt',
                      name:'units',
                      data:'cm/s'}

    plasma_vz_desc = {attribute,obj:'/plasma/vz',
                      name:'description',
                      data:'Bulk plasma flow in the z-direction: Vz(r,z)'}
    plasma_vz_unit = {attribute,obj:'/plasma/vz',
                      name:'units',
                      data:'cm/s'}

    plasma_nr_desc = {attribute,obj:'/plasma/nr',
                      name:'description',
                      data:'Number of R values'}

    plasma_nz_desc = {attribute,obj:'/plasma/nz',
                      name:'description',
                      data:'Number of Z values'}

    plasma_r_desc = {attribute,obj:'/plasma/r',
                     name:'description',
                     data:'Radius'}
    plasma_r_unit = {attribute,obj:'/plasma/r',
                     name:'units',
                     data:'cm'}

    plasma_z_desc = {attribute,obj:'/plasma/z',
                     name:'description',
                     data:'Z'}
    plasma_z_unit = {attribute,obj:'/plasma/z',
                     name:'units',
                     data:'cm'}

    plasma_r2d_desc = {attribute,obj:'/plasma/r2d',
                       name:'description',
                       data:'Radius grid: R(r,z)'}
    plasma_r2d_unit = {attribute,obj:'/plasma/r2d',
                       name:'units',
                       data:'cm'}

    plasma_z2d_desc = {attribute,obj:'/plasma/z2d',
                       name:'description',
                       data:'Z grid: Z(r,z)'}
    plasma_z2d_unit = {attribute,obj:'/plasma/z2d',
                       name:'units',
                       data:'cm'}

    plasma_mask_desc = {attribute,obj:'/plasma/mask',
                         name:'description',
                         data:'Boolean mask that indicates where' +
                         ' the plasma parameters are well defined'}

    plasma_atts = [plasma_desc, plasma_cs, plasma_ds_desc,
                   plasma_time_desc, plasma_time_unit,
                   plasma_mask_desc,
                   plasma_dene_desc, plasma_dene_unit,
                   plasma_te_desc, plasma_te_unit,
                   plasma_ti_desc, plasma_ti_unit,
                   plasma_zeff_desc,
                   plasma_vr_desc, plasma_vr_unit,
                   plasma_vt_desc, plasma_vt_unit,
                   plasma_vz_desc, plasma_vz_unit,
                   plasma_nr_desc, plasma_nz_desc,
                   plasma_r_desc, plasma_r_unit,
                   plasma_z_desc, plasma_z_unit,
                   plasma_r2d_desc, plasma_r2d_unit,
                   plasma_z2d_desc, plasma_z2d_unit ]

    # Electromagnetic fields attributes
    fields_desc = {attribute,obj:'/fields',
                   name:'description',
                   data:'Electromagnetic Fields'}

    fields_cs = {attribute,obj:'/fields',
                 name:'coordinate_system',
                 data:'Cylindrical'}

    fields_ds_desc = {attribute,obj:'/fields/data_source',
                      name:'description',
                      data:'Source of the EM equilibrium'}

    fields_mask_desc = {attribute,obj:'/fields/mask',
                        name:'description',
                        data:'Boolean mask that indicates where' +
                        ' the fields are well defined'}

    fields_time_desc = {attribute,obj:'/fields/time',
                        name:'description',
                        data:'Time'}
    fields_time_unit = {attribute,obj:'/fields/time',
                        name:'units',
                        data:'s'}

    fields_br_desc = {attribute,obj:'/fields/br',
                      name:'description',
                      data:'Magnetic field in the r-direction: Br(r,z)'}
    fields_br_unit = {attribute,obj:'/fields/br',
                      name:'units',
                      data:'T'}

    fields_bt_desc = {attribute,obj:'/fields/bt',
                      name:'description',
                      data:'Magnetic field in the theta/torodial-direction: Bt(r,z)'}
    fields_bt_unit = {attribute,obj:'/fields/bt',
                      name:'units',
                      data:'T'}

    fields_bz_desc = {attribute,obj:'/fields/bz',
                      name:'description',
                      data:'Magnetic field in the z-direction: Bz(r,z)'}
    fields_bz_unit = {attribute,obj:'/fields/bz',
                      name:'units',
                      data:'T'}

    fields_er_desc = {attribute,obj:'/fields/er',
                      name:'description',
                      data:'Electric field in the r-direction: Er(r,z)'}
    fields_er_unit = {attribute,obj:'/fields/er',
                      name:'units',
                      data:'V/m'}

    fields_et_desc = {attribute,obj:'/fields/et',
                      name:'description',
                      data:'Electric field in the theta/torodial-direction: Et(r,z)'}
    fields_et_unit = {attribute,obj:'/fields/et',
                      name:'units',
                      data:'V/m'}

    fields_ez_desc = {attribute,obj:'/fields/ez',
                      name:'description',
                      data:'Electric field in the z-direction: Ez(r,z)'}
    fields_ez_unit = {attribute,obj:'/fields/ez',
                      name:'units',
                      data:'V/m'}

    fields_nr_desc = {attribute,obj:'/fields/nr',
                      name:'description',
                      data:'Number of R values'}

    fields_nz_desc = {attribute,obj:'/fields/nz',
                      name:'description',
                      data:'Number of Z values'}

    fields_r_desc = {attribute,obj:'/fields/r',
                     name:'description',
                     data:'Radius'}
    fields_r_unit = {attribute,obj:'/fields/r',
                     name:'units',
                     data:'cm'}

    fields_z_desc = {attribute,obj:'/fields/z',
                     name:'description',
                     data:'Z'}
    fields_z_unit = {attribute,obj:'/fields/z',
                     name:'units',
                     data:'cm'}

    fields_r2d_desc = {attribute,obj:'/fields/r2d',
                       name:'description',
                       data:'Radius grid: R(r,z)'}
    fields_r2d_unit = {attribute,obj:'/fields/r2d',
                       name:'units',
                       data:'cm'}

    fields_z2d_desc = {attribute,obj:'/fields/z2d',
                       name:'description',
                       data:'Z grid: Z(r,z)'}
    fields_z2d_unit = {attribute,obj:'/fields/z2d',
                       name:'units',
                       data:'cm'}

    fields_atts = [fields_desc, fields_cs, fields_ds_desc,
                   fields_mask_desc,
                   fields_time_desc, fields_time_unit,
                   fields_br_desc, fields_br_unit,
                   fields_bt_desc, fields_bt_unit,
                   fields_bz_desc, fields_bz_unit,
                   fields_er_desc, fields_er_unit,
                   fields_et_desc, fields_et_unit,
                   fields_ez_desc, fields_ez_unit,
                   fields_nr_desc, fields_nz_desc,
                   fields_r_desc, fields_r_unit,
                   fields_z_desc, fields_z_unit,
                   fields_r2d_desc, fields_r2d_unit,
                   fields_z2d_desc, fields_z2d_unit ]

    atts = [root_atts, plasma_atts, fields_atts]

    write_hdf5,["plasma","fields"],filename=filename, atts=atts, /clobber

    if file_test(filename):
        success('Equilibrium file created: '+filename
     else:
        error('Equilibrium file creation failed.'



def write_distribution, filename, distri):
    """
    ##`write_distribution, filename, dist`
    Write fast-ion distribution to a HDF5 file

    ###Input Arguments
         **filename**: Name of the distribution file

         **dist**: Fast-ion distribution structure

    ###Example Usage
    ```idl
    IDL> write_distribution, filename, distri
    ```
    """
    info('Writing fast-ion distribution file...'

    root_atts = {attribute,obj:'/',
                 name:'description',
                 data:'Fast-ion distribution for FIDASIM'}

    cs_desc = {attribute,obj:'/',
               name:'coordinate_system',
               data:'Cylindrical'}

    ds_desc = {attribute,obj:'/data_source',
               name:'description',
               data:'Source of the fast-ion distribution'}

    type_desc = {attribute,obj:'/type',
                 name:'description',
                 data:'Distribution type: '+
                      '1="Guiding Center Density Function", '+
                      '2="Guiding Center Monte Carlo", '+
                      '3="Full Orbit Monte Carlo"'}

    time_desc = {attribute,obj:'/time',
                 name:'description',
                 data:'Distribution time'}
    time_unit = {attribute,obj:'/time',
                 name:'units',
                 data:'s'}

    type = distri.type

    if type eq 1:
        nen_desc = {attribute,obj:'/nenergy',
                    name:'description',
                    data:'Number of energy values'}
        np_desc = {attribute,obj:'/npitch',
                   name:'description',
                   data:'Number of pitch values'}

        energy_desc = {attribute,obj:'/energy',
                       name:'description',
                       data:'Energy'}
        energy_unit = {attribute,obj:'/energy',
                       name:'units',
                       data:'keV'}

        pitch_desc = {attribute,obj:'/pitch',
                       name:'description',
                       data:'Pitch: p = v_parallel/v  w.r.t. the magnetic field'}

        f_desc = {attribute,obj:'/f',
                  name:'description',
                  data:'Fast-ion density function: F(E,p,R,Z)'}
        f_unit = {attribute,obj:'/f',
                  name:'units',
                  data:'fast-ions/(dE*dP*cm^3)'}

        denf_desc = {attribute,obj:'/denf',
                     name:'description',
                     data:'Fast-ion density: Denf(r,z)'}
        denf_unit = {attribute,obj:'/denf',
                     name:'units',
                     data:'cm^-3'}

        nr_desc = {attribute,obj:'/nr',
                   name:'description',
                   data:'Number of R values'}

        nz_desc = {attribute,obj:'/nz',
                   name:'description',
                   data:'Number of Z values'}

        r_desc = {attribute,obj:'/r',
                  name:'description',
                  data:'Radius'}
        r_unit = {attribute,obj:'/r',
                  name:'units',
                  data:'cm'}

        z_desc = {attribute,obj:'/z',
                  name:'description',
                  data:'Z'}
        z_unit = {attribute,obj:'/z',
                  name:'units',
                  data:'cm'}

        r2d_desc = {attribute,obj:'/r2d',
                    name:'description',
                    data:'Radius grid: R(r,z)'}
        r2d_unit = {attribute,obj:'/r2d',
                    name:'units',
                    data:'cm'}

        z2d_desc = {attribute,obj:'/z2d',
                    name:'description',
                    data:'Z grid: Z(r,z)'}
        z2d_unit = {attribute,obj:'/z2d',
                    name:'units',
                    data:'cm'}

        atts = [root_atts, cs_desc, ds_desc,
                type_desc, time_desc,time_unit,
                nen_desc, np_desc,
                energy_desc, energy_unit,
                pitch_desc,
                f_desc, f_unit,
                denf_desc, denf_unit,
                nr_desc, nz_desc,
                r_desc, r_unit,
                z_desc, z_unit,
                r2d_desc, r2d_unit,
                z2d_desc, z2d_unit]

     else:

        np_desc = {attribute,obj:'/nparticle',
                   name:'description',
                   data:'Number of MC particles'}

        nc_desc = {attribute,obj:'/nclass',
                   name:'description',
                   data:'Number of orbit classes'}

        r_desc = {attribute,obj:'/r',
                  name:'description',
                  data:'R position of a MC particle'}
        r_unit = {attribute,obj:'/r',
                  name:'units',
                  data:'cm'}

        z_desc = {attribute,obj:'/z',
                  name:'description',
                  data:'Z position of a MC particle'}
        z_unit = {attribute,obj:'/z',
                  name:'units',
                  data:'cm'}

        w_desc = {attribute,obj:'/weight',
                  name:'description',
                  data:'Weight of a MC particle: sum(weight) = # of fast-ions '}
        w_unit = {attribute,obj:'/weight',
                  name:'units',
                  data:'fast-ions/particle'}

        c_desc = {attribute,obj:'/class',
                  name:'description',
                  data:'Orbit class of a MC particle: class in Set(1:nclass)'}

        if type eq 2:
            energy_desc = {attribute,obj:'/energy',
                           name:'description',
                           data:'Energy of a MC particle'}
            energy_unit = {attribute,obj:'/energy',
                           name:'units',
                           data:'keV'}

            pitch_desc = {attribute,obj:'/pitch',
                          name:'description',
                          data:'Pitch of a MC particle: p = v_parallel/v  w.r.t. the magnetic field'}
            type_atts = [energy_desc,energy_unit,pitch_desc]
         else:
            vr_desc = {attribute,obj:'/vr',
                       name:'description',
                       data:'Radial velocity of a MC particle'}
            vr_unit = {attribute,obj:'/vr',
                       name:'units',
                       data:'cm/s'}
            vt_desc = {attribute,obj:'/vt',
                       name:'description',
                       data:'Torodial velocity of a MC particle'}
            vt_unit = {attribute,obj:'/vt',
                       name:'units',
                       data:'cm/s'}
            vz_desc = {attribute,obj:'/vz',
                       name:'description',
                       data:'Z velocity of a MC particle'}
            vz_unit = {attribute,obj:'/vz',
                       name:'units',
                       data:'cm/s'}
            type_atts = [vr_desc,vr_unit,vt_desc,vt_unit,vz_desc,vz_unit]


        atts = [root_atts, cs_desc, ds_desc,
                type_desc, time_desc,time_unit,
                np_desc, nc_desc,
                r_desc, r_unit,
                z_desc, z_unit,
                w_desc, w_unit,
                c_desc, type_atts]


    write_hdf5, distri, filename=filename,atts=atts, /clobber, compress=4

    if file_test(filename):
        success('Distribution file created: ' + filename)
    else:
        error('Distribution file creation failed.')



