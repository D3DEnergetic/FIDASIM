PRO check_inputs, inputs, err_status
    ;+##`check_inputs, inputs, err`
    ;+Checks if input structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: input structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> check_inputs, inputs, err
    ;+```
    info,'Checking simulation settings...'
    err_status = 0
    zero_string = {dims:0,type:'STRING'}
    zero_int = {dims:0,type:'INT'}
    zero_long = {dims:0,type:'LONG'}
    zero_double = {dims:0,type:'DOUBLE'}
    three_double = {dims:[3],type:'DOUBLE'}
    schema = {comment:zero_string, $
              shot:zero_long, time:zero_double, $
              runid:zero_string, device:zero_string, $ 
              install_dir:zero_string, tables_file:zero_string, result_dir:zero_string, $
              nlambda:zero_int, lambdamin:zero_double, lambdamax:zero_double, $
              nx:zero_int, ny:zero_int, nz:zero_int, $
              alpha:zero_double, beta:zero_double, gamma:zero_double, $
              origin:three_double, xmin:zero_double, xmax:zero_double, $
              ymin:zero_double, ymax:zero_double, zmin:zero_double, zmax:zero_double, $
              ab:zero_double, ai:zero_double, species_mix:three_double, $
              pinj:zero_double, einj:zero_double, impurity_charge:zero_int, $
              n_fida:zero_long, n_nbi:zero_long, n_dcx:zero_long, $
              n_npa:zero_long, n_halo:zero_long, n_birth:zero_long, $
              ne_wght:zero_int, np_wght:zero_int, nphi_wght:zero_int, $
              emax_wght:zero_double, nlambda_wght:zero_int, $
              lambdamin_wght:zero_double, lambdamax_wght:zero_double, $
              calc_npa:zero_int, calc_fida:zero_int, calc_bes:zero_int, $
              calc_brems:zero_int, calc_birth:zero_int, $
              calc_fida_wght:zero_int, calc_npa_wght:zero_int, $
              dump_dcx:zero_int, load_neutrals:zero_int, verbose:zero_int}

    check_struct_schema, schema, inputs, err_status, desc="simulation settings"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    ;Normalize File Paths
    inputs.result_dir = expand_path(inputs.result_dir)
    inputs.install_dir = expand_path(inputs.install_dir)

    if inputs.alpha gt 2*!DPI or $
       inputs.beta gt 2*!DPI or $
       inputs.gamma gt 2*!DPI then begin
        error,'Angles must be in radians'
        err_status = 1
    endif

    if inputs.lambdamin ge inputs.lambdamax then begin
        error,'Invalid wavelength range. Expected lambdamin < lamdbdamax'
        err_status = 1
    endif

    if inputs.lambdamin_wght ge inputs.lambdamax_wght then begin
        error,'Invalid wavelength range. Expected lambdamin_wght < lamdbdamax_wght'
        err_status = 1
    endif

    if inputs.xmin ge inputs.xmax then begin
        error,'Invalid x range. Expected xmin < xmax'
        err_status = 1
    endif

    if inputs.ymin ge inputs.ymax then begin
        error,'Invalid y range. Expected ymin < ymax'
        err_status = 1
    endif

    if inputs.zmin ge inputs.zmax then begin
        error,'Invalid z range. Expected zmin < zmax'
        err_status = 1
    endif

    if inputs.pinj le 0. or inputs.einj le 0.0 then begin
        error,'The selected source is not on'
        print,'einj = ',inputs.einj
        print,'pinj = ',inputs.pinj
        err_status = 1
    endif

    if abs(total(inputs.species_mix) - 1.0) gt 1.d-3 then begin
        error,'species_mix does not sum to 1.0'
        print,'sum(species_mix) = ',total(inputs.species_mix)
        err_status = 1
    endif

    if inputs.impurity_charge le 1 then begin
        error,'Invalid impurity charge. Expected impurity charge > 1'
        err_status = 1
    endif

    ps = path_sep()
    input_file = inputs.result_dir+ps+inputs.runid+'_inputs.dat'
    equilibrium_file = inputs.result_dir+ps+inputs.runid+'_equilibrium.h5'
    geometry_file = inputs.result_dir+ps+inputs.runid+'_geometry.h5'
    distribution_file = inputs.result_dir+ps+inputs.runid+'_distribution.h5'
    neutrals_file = inputs.result_dir+ps+inputs.runid+'_neutrals.h5'

    inputs = create_struct(inputs,'input_file',input_file, $
                                  'equilibrium_file',equilibrium_file,$ 
                                  'geometry_file',geometry_file, $
                                  'distribution_file',distribution_file, $
                                  'neutrals_file',neutrals_file)

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid simulation settings. Exiting...'
    endif else begin
        success,'Simulation settings are valid'
    endelse

END
