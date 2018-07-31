PRO write_distribution, filename, distri
    ;+##`write_distribution, filename, dist`
    ;+Write fast-ion distribution to a HDF5 file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the distribution file
    ;+
    ;+     **dist**: Fast-ion distribution structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_distribution, filename, distri
    ;+```

    info, 'Writing fast-ion distribution file...'

    root_atts = {attribute,obj:'/', $
                 name:'description', $
                 data:'Fast-ion distribution for FIDASIM'}

    cs_desc = {attribute,obj:'/', $
               name:'coordinate_system', $
               data:'Cylindrical'}

    ds_desc = {attribute,obj:'/data_source', $
               name:'description', $
               data:'Source of the fast-ion distribution'}    
    
    type_desc = {attribute,obj:'/type', $
                 name:'description', $
                 data:'Distribution type: '+ $
                      '1="Guiding Center Density Function", '+ $
                      '2="Guiding Center Monte Carlo", '+ $
                      '3="Full Orbit Monte Carlo"'}

    time_desc = {attribute,obj:'/time', $
                 name:'description', $
                 data:'Distribution time'}
    time_unit = {attribute,obj:'/time', $
                 name:'units', $
                 data:'s'}

    type = distri.type

    if type eq 1 then begin
        nen_desc = {attribute,obj:'/nenergy',$
                    name:'description', $
                    data:'Number of energy values'}
        np_desc = {attribute,obj:'/npitch', $
                   name:'description', $
                   data:'Number of pitch values'}

        energy_desc = {attribute,obj:'/energy', $
                       name:'description',$
                       data:'Energy'}
        energy_unit = {attribute,obj:'/energy', $
                       name:'units', $
                       data:'keV'}
 
        pitch_desc = {attribute,obj:'/pitch', $
                       name:'description',$
                       data:'Pitch: p = v_parallel/v  w.r.t. the magnetic field'}

        f_desc = {attribute,obj:'/f', $
                  name:'description', $
                  data:'Fast-ion density function: F(E,p,R,Z,Phi)'}
        f_unit = {attribute,obj:'/f', $
                  name:'units', $
                  data:'fast-ions/(dE*dP*cm^3)'} 

        denf_desc = {attribute,obj:'/denf', $
                     name:'description', $
                     data:'Fast-ion density: Denf(r,z,phi)'}
        denf_unit = {attribute,obj:'/denf', $
                     name:'units', $
                     data:'cm^-3'} 

        nr_desc = {attribute,obj:'/nr', $
                   name:'description', $
                   data:'Number of R values'}
   
        nz_desc = {attribute,obj:'/nz', $
                   name:'description', $
                   data:'Number of Z values'}

        nphi_desc = {attribute,obj:'/nphi', $
                   name:'description', $
                   data:'Number of Phi values'}

        r_desc = {attribute,obj:'/r', $
                  name:'description', $
                  data:'Radius'}
        r_unit = {attribute,obj:'/r', $
                  name:'units', $
                  data:'cm'}

        z_desc = {attribute,obj:'/z', $
                  name:'description', $
                  data:'Z'}
        z_unit = {attribute,obj:'/z', $
                  name:'units', $
                  data:'cm'}

        phi_desc = {attribute,obj:'/phi', $
                  name:'description', $
                  data:'Phi'}
        phi_unit = {attribute,obj:'/phi', $
                  name:'units', $
                  data:'rad'}

        r2d_desc = {attribute,obj:'/r2d', $
                    name:'description', $
                    data:'Radius grid: R(r,z)'}
        r2d_unit = {attribute,obj:'/r2d', $
                    name:'units', $
                    data:'cm'}

        z2d_desc = {attribute,obj:'/z2d', $
                    name:'description', $
                    data:'Z grid: Z(r,z)'}
        z2d_unit = {attribute,obj:'/z2d', $
                    name:'units', $
                    data:'cm'}

        atts = [root_atts, cs_desc, ds_desc, $
                type_desc, time_desc,time_unit, $
                nen_desc, np_desc, $
                energy_desc, energy_unit, $
                pitch_desc, $
                f_desc, f_unit, $
                denf_desc, denf_unit, $
                nr_desc, nz_desc, $
                nphi_desc, $
                r_desc, r_unit, $
                z_desc, z_unit, $
                phi_desc, phi_unit, $
                r2d_desc, r2d_unit, $
                z2d_desc, z2d_unit]
        
    endif else begin

        np_desc = {attribute,obj:'/nparticle', $
                   name:'description', $
                   data:'Number of MC particles'}

        nc_desc = {attribute,obj:'/nclass', $
                   name:'description', $
                   data:'Number of orbit classes'}

        r_desc = {attribute,obj:'/r', $
                  name:'description', $
                  data:'R position of a MC particle'}
        r_unit = {attribute,obj:'/r', $
                  name:'units', $
                  data:'cm'}

        z_desc = {attribute,obj:'/z', $
                  name:'description', $
                  data:'Z position of a MC particle'}
        z_unit = {attribute,obj:'/z', $
                  name:'units', $
                  data:'cm'}

        phi_desc = {attribute,obj:'/phi', $
                  name:'description', $
                  data:'Phi position of a MC particle'}
        phi_unit = {attribute,obj:'/phi', $
                  name:'units', $
                  data:'rad'}

        w_desc = {attribute,obj:'/weight', $
                  name:'description', $
                  data:'Weight of a MC particle: sum(weight) = # of fast-ions '}
        w_unit = {attribute,obj:'/weight', $
                  name:'units', $
                  data:'fast-ions/particle'}

        c_desc = {attribute,obj:'/class', $
                  name:'description', $
                  data:'Orbit class of a MC particle: class in Set(1:nclass)'}

        if type eq 2 then begin
            energy_desc = {attribute,obj:'/energy', $
                           name:'description', $
                           data:'Energy of a MC particle'}
            energy_unit = {attribute,obj:'/energy', $
                           name:'units', $
                           data:'keV'}

            pitch_desc = {attribute,obj:'/pitch', $
                          name:'description', $
                          data:'Pitch of a MC particle: p = v_parallel/v  w.r.t. the magnetic field'}
            type_atts = [energy_desc,energy_unit,pitch_desc]
        endif else begin
            vr_desc = {attribute,obj:'/vr', $
                       name:'description', $
                       data:'Radial velocity of a MC particle'}
            vr_unit = {attribute,obj:'/vr', $
                       name:'units', $
                       data:'cm/s'}
            vt_desc = {attribute,obj:'/vt', $
                       name:'description', $
                       data:'Torodial velocity of a MC particle'}
            vt_unit = {attribute,obj:'/vt', $
                       name:'units', $
                       data:'cm/s'}
            vz_desc = {attribute,obj:'/vz', $
                       name:'description', $
                       data:'Z velocity of a MC particle'}
            vz_unit = {attribute,obj:'/vz', $
                       name:'units', $
                       data:'cm/s'}
            type_atts = [vr_desc,vr_unit,vt_desc,vt_unit,vz_desc,vz_unit]
        endelse

        atts = [root_atts, cs_desc, ds_desc, $
                type_desc, time_desc,time_unit, $
                np_desc, nc_desc, $
                r_desc, r_unit, $
                z_desc, z_unit, $
                phi_desc, phi_unit, $
                w_desc, w_unit, $
                c_desc, type_atts] 
    endelse

    write_hdf5, distri, filename=filename,atts=atts, /clobber, compress=4
   
    if file_test(filename) then begin
        success, 'Distribution file created: '+filename
    endif else begin
        error, 'Distribution file creation failed.'
    endelse
END

