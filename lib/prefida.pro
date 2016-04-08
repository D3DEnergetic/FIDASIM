PRO prefida,inputs,grid,nbi,plasma,fields,fbm,spec=spec,npa=npa
    ;+##`prefida, inputs, grid, nbi, plasma, fields, dist, spec=spec, npa=npa`
    ;+Checks FIDASIM inputs and writes FIDASIM input files
    ;+
    ;+###Input Arguments
    ;+     **inputs**: Inputs structure
    ;+
    ;+     **grid**: Interpolation grid structure
    ;+
    ;+     **nbi**: Neutral beam geometry structure
    ;+
    ;+     **plasma**: Plasma parameters structure
    ;+
    ;+     **fields**: Electromagnetic fields structure
    ;+
    ;+     **dist**: Fast-ion distribution structure
    ;+
    ;+###Keyword Arguments
    ;+     **spec**: Optional, Spectral geometry structure
    ;+
    ;+     **npa**: Optional, NPA geometry structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> prefida, inputs, grid, nbi, plasma, fields, dist, spec=spec, npa=npa
    ;+```
    COMPILE_OPT DEFINT32

    ;;CHECK INPUTS
    check_inputs,inputs,err_status
    if err_status ne 0 then goto, GET_OUT

    ;;MAKE DIRECTORIES IF THEY DONT EXIST
    if file_test(inputs.result_dir,/directory) eq 0 then begin
        file_mkdir,inputs.result_dir
    endif

    ;;CHECK INTERPOLATION GRID
    check_grid, grid, err_status
    if err_status ne 0 then goto, GET_OUT

    ;;CHECK BEAM INPUTS
    check_beam, inputs, nbi, err_status
    if err_status ne 0 then goto, GET_OUT

    ;;CHECK PLASMA PARAMETERS
    check_plasma, inputs, grid, plasma, err_status
    if err_status ne 0 then goto, GET_OUT

    ;;CHECK ELECTROMAGNETIC FIELDS
    check_fields, inputs, grid, fields, err_status
    if err_status ne 0 then goto, GET_OUT

    ;;CHECK FAST-ION DISTRIBUTION
    check_distribution, inputs, grid, fbm, err_status
    if err_status ne 0 then goto, GET_OUT

    ;;CHECK FIDA/BES
    if keyword_set(spec) then begin
        check_spec, inputs, spec, err_status
        if err_status ne 0 then goto, GET_OUT
    endif

    ;;CHECK NPA
    if keyword_set(npa) then begin
        check_npa, inputs, npa, err_status
        if err_status ne 0 then goto, GET_OUT
    endif

    ;;WRITE FIDASIM INPUT FILES
    write_namelist, inputs.input_file, inputs

    ;;WRITE GEOMETRY FILE
    write_geometry, inputs.geometry_file, nbi, spec=spec, npa=npa

    ;;WRITE EQUILIBRIUM FILE
    write_equilibrium, inputs.equilibrium_file, plasma, fields
  
    ;;WRITE DISTRIBUTION FILE
    write_distribution, inputs.distribution_file, fbm

    print,''
    print,''
    success,'FIDASIM pre-processing completed'
    print, 'To run FIDASIM use the following command'
    print, inputs.install_dir+'/fidasim '+inputs.result_dir+'/'+inputs.runid+'_inputs.dat'
    print,''
    print,''
    GET_OUT:
END
