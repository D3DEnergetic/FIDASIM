PRO check_fields, inp, grid, fields
    ;+#check_fields
    ;+Checks if electromagnetic fields structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: Input structure
    ;+ 
    ;+     **grid**: Interpolation grid structure
    ;+ 
    ;+     **fields**: Electromagnetic fields structure
    ;+ 
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_fields, inputs, grid, fields
    ;+```

    err_status=0
    info,'Checking electromagnetic fields...'

    nr = grid.nr
    nz = grid.nz
    w = where("nphi" eq strlowcase(TAG_NAMES(grid)),nw)
    if nw eq 0 then begin
        info,'"nphi" is missing from the fields, assuming axisymmetry'
        nphi = 1
    endif else begin
        nphi = grid.nphi
    endelse

    zero_string = {dims:0, type:'STRING'}
    zero_double = {dims:0, type:'DOUBLE'}
    nrnznphi_double = {dims:[nr,nz,nphi], type:'DOUBLE'}
    nrnznphi_int = {dims:[nr,nz,nphi], type:'INT'}
    schema = {time:zero_double, $
              br:nrnznphi_double, $
              bt:nrnznphi_double, $
              bz:nrnznphi_double, $
              er:nrnznphi_double, $
              et:nrnznphi_double, $
              ez:nrnznphi_double, $
              mask:nrnznphi_int, $
              data_source:zero_string}

    check_struct_schema,schema,fields,err_status, desc = "electromagnetic fields"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    if fields.data_source eq '' then begin
        error, 'Invalid data source. An empty string is not a data source.'
        err_status = 1
    endif
 
    if abs(fields.time - inp.time) gt 0.02 then begin
        warn,'Electromagnetic fields time and input time do not match'
        print,'Input time: ',inp.time
        print,'Electromagnetic fields time: ',fields.time
    endif

    fields = create_struct(fields, grid)
    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid electromagnetic fields. Exiting...',/halt
    endif else begin
        success,'Electromagnetic fields are valid'
    endelse
END
