PRO check_inputs, inputs
    ;+#check_inputs
    ;+Checks if input structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+ 
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_inputs, inputs
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
              tables_file:zero_string, result_dir:zero_string, $
              nlambda:zero_int, lambdamin:zero_double, lambdamax:zero_double, $
              nx:zero_int, ny:zero_int, nz:zero_int, $
              alpha:zero_double, beta:zero_double, gamma:zero_double, $
              origin:three_double, xmin:zero_double, xmax:zero_double, $
              ymin:zero_double, ymax:zero_double, zmin:zero_double, zmax:zero_double, $
              ab:zero_double, current_fractions:three_double, $
              pinj:zero_double, einj:zero_double, $
              n_fida:zero_long, n_pfida:zero_long, n_nbi:zero_long, n_dcx:zero_long, $
              n_npa:zero_long, n_pnpa:zero_long,n_halo:zero_long, n_birth:zero_long, $
              ne_wght:zero_int, np_wght:zero_int, nphi_wght:zero_int, $
              emax_wght:zero_double, nlambda_wght:zero_int, $
              lambdamin_wght:zero_double, lambdamax_wght:zero_double, $
              calc_npa:zero_int, calc_fida:zero_int, calc_bes:zero_int, $
              calc_dcx:zero_int, calc_halo:zero_int, calc_cold:zero_int, $
	      calc_pnpa:zero_int, calc_pfida:zero_int, $
              calc_brems:zero_int, calc_birth:zero_int, calc_neutron:zero_int,$
              calc_fida_wght:zero_int, calc_npa_wght:zero_int}

    check_struct_schema, schema, inputs, err_status, desc="simulation settings"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    ;Normalize File Paths
    inputs.result_dir = expand_path(inputs.result_dir)

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

    if abs(total(inputs.current_fractions) - 1.0) gt 1.d-3 then begin
        error,'current_fractions do not sum to 1.0'
        print,'sum(current_fractions) = ',total(inputs.current_fractions)
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
                                  'load_neutrals', 0, $
                                  'flr',2, $
                                  'verbose', 1, $
                                  'seed', -1, $
                                  'stark_components', 0, $
                                  'neutrals_file',neutrals_file)
    if (total(tag_names(inputs) EQ 'OUTPUT_NEUTRAL_RESERVOIR') EQ 0) then inputs=create_struct(inputs,'output_neutral_reservoir',1)

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid simulation settings. Exiting...',/halt
    endif else begin
        success,'Simulation settings are valid'
    endelse

END
