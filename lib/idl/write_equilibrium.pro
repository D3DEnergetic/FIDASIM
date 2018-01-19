PRO write_equilibrium, filename, plasma, fields
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

    info, 'Writing equilibrium file...'

    root_atts = {attribute,obj:'/', $
                 name:'description', $
                 data:'Plasma Parameters and Electromagnetic Fields for FIDASIM'}

    ;; Plasma Attributes
    plasma_desc = {attribute,obj:'/plasma', $
                   name:'description', $
                   data:'Plasma Parameters'}

    plasma_cs = {attribute,obj:'/plasma', $
                 name:'coordinate_system', $
                 data:'Cylindrical'}

    plasma_ds_desc = {attribute, obj:'/plasma/data_source', $
                      name:'description', $
                      data:'Source of the plasma parameters'}

    plasma_time_desc = {attribute,obj:'/plasma/time', $
                        name:'description', $
                        data:'Time'}
    plasma_time_unit = {attribute,obj:'/plasma/time', $
                        name:'units', $
                        data:'s'}

    plasma_dene_desc = {attribute,obj:'/plasma/dene', $
                        name:'description', $
                        data:'Electron Number Density: Dene(r,z)'}
    plasma_dene_unit = {attribute,obj:'/plasma/dene', $
                        name:'units', $
                        data:'cm^-3'}

    plasma_denn_desc = {attribute,obj:'/plasma/denn', $
                        name:'description', $
                        data:'Cold Neutral Number Density: Denn(r,z)'}
    plasma_denn_unit = {attribute,obj:'/plasma/denn', $
                        name:'units', $
                        data:'cm^-3'}

    plasma_te_desc = {attribute,obj:'/plasma/te', $
                      name:'description', $
                      data:'Electron Temperature: Te(r,z)'}
    plasma_te_unit = {attribute,obj:'/plasma/te', $
                      name:'units', $
                      data:'keV'}

    plasma_ti_desc = {attribute,obj:'/plasma/ti', $
                      name:'description', $
                      data:'Ion Temperature: Ti(r,z)'}
    plasma_ti_unit = {attribute,obj:'/plasma/ti', $
                      name:'units', $
                      data:'keV'}

    plasma_zeff_desc = {attribute,obj:'/plasma/zeff', $
                        name:'description', $
                        data:'Effective Nuclear Charge: Zeff(r,z)'}

    plasma_vr_desc = {attribute,obj:'/plasma/vr', $
                      name:'description', $
                      data:'Bulk plasma flow in the r-direction: Vr(r,z)'}
    plasma_vr_unit = {attribute,obj:'/plasma/vr', $
                      name:'units', $
                      data:'cm/s'}

    plasma_vt_desc = {attribute,obj:'/plasma/vt', $
                      name:'description', $
                      data:'Bulk plasma flow in the theta/torodial-direction: Vt(r,z)'}
    plasma_vt_unit = {attribute,obj:'/plasma/vt', $
                      name:'units', $
                      data:'cm/s'}

    plasma_vz_desc = {attribute,obj:'/plasma/vz', $
                      name:'description', $
                      data:'Bulk plasma flow in the z-direction: Vz(r,z)'}
    plasma_vz_unit = {attribute,obj:'/plasma/vz', $
                      name:'units', $
                      data:'cm/s'}

    plasma_nr_desc = {attribute,obj:'/plasma/nr', $
                      name:'description', $
                      data:'Number of R values'}

    plasma_nz_desc = {attribute,obj:'/plasma/nz', $
                      name:'description', $
                      data:'Number of Z values'}

    plasma_r_desc = {attribute,obj:'/plasma/r', $
                     name:'description', $
                     data:'Radius'}
    plasma_r_unit = {attribute,obj:'/plasma/r', $
                     name:'units', $
                     data:'cm'}

    plasma_z_desc = {attribute,obj:'/plasma/z', $
                     name:'description', $
                     data:'Z'}
    plasma_z_unit = {attribute,obj:'/plasma/z', $
                     name:'units', $
                     data:'cm'}

    plasma_r2d_desc = {attribute,obj:'/plasma/r2d', $
                       name:'description', $
                       data:'Radius grid: R(r,z)'}
    plasma_r2d_unit = {attribute,obj:'/plasma/r2d', $
                       name:'units', $
                       data:'cm'}

    plasma_z2d_desc = {attribute,obj:'/plasma/z2d', $
                       name:'description', $
                       data:'Z grid: Z(r,z)'}
    plasma_z2d_unit = {attribute,obj:'/plasma/z2d', $
                       name:'units', $
                       data:'cm'}

    plasma_mask_desc = {attribute,obj:'/plasma/mask',$
                         name:'description', $
                         data:'Boolean mask that indicates where' +$ 
                         ' the plasma parameters are well defined'}

    plasma_atts = [plasma_desc, plasma_cs, plasma_ds_desc, $
                   plasma_time_desc, plasma_time_unit, $
                   plasma_mask_desc, $
                   plasma_dene_desc, plasma_dene_unit, $
                   plasma_denn_desc, plasma_denn_unit, $
                   plasma_te_desc, plasma_te_unit, $
                   plasma_ti_desc, plasma_ti_unit, $
                   plasma_zeff_desc, $
                   plasma_vr_desc, plasma_vr_unit, $
                   plasma_vt_desc, plasma_vt_unit, $
                   plasma_vz_desc, plasma_vz_unit, $
                   plasma_nr_desc, plasma_nz_desc, $
                   plasma_r_desc, plasma_r_unit, $
                   plasma_z_desc, plasma_z_unit, $
                   plasma_r2d_desc, plasma_r2d_unit, $
                   plasma_z2d_desc, plasma_z2d_unit ]

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
