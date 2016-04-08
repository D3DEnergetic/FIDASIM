PRO check_plasma, inp, grid, plasma, err_status
    ;+#check_plasma
    ;+Checks if plasma paramters structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: Input structure
    ;+
    ;+     **grid**: Interpolation grid structure
    ;+
    ;+     **plasma**: Plasma parameters structure
    ;+ 
    ;+##Output Arguments
    ;+     **err**: error code
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_plasma, inputs, grid, plasma, err
    ;+```
    err_status=0
    info,'Checking plasma parameters...'

    nr = grid.nr
    nz = grid.nz
    zero_string = {dims:0, type:'STRING'}
    zero_double = {dims:0, type:'DOUBLE'}
    nrnz_double = {dims:[nr,nz], type:'DOUBLE'}
    nrnz_int    = {dims:[nr,nz], type:'INT'}
    schema = {time:zero_double, $
              vr:nrnz_double, $
              vt:nrnz_double, $
              vz:nrnz_double, $
              dene:nrnz_double, $
              ti:nrnz_double, $
              te:nrnz_double, $
              zeff:nrnz_double, $
              mask:nrnz_int, $
              data_source:zero_string}

    check_struct_schema,schema,plasma,err_status,desc = "plasma parameters"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif
  
    if plasma.data_source eq '' then begin
        error, 'Invalid data source. An empty string is not a data source.'
        err_status = 1
    endif
 
    ;;Electron density
    plasma.dene = plasma.dene > 0. ;[1/cm^3]

    ;;Zeff
    plasma.zeff = plasma.zeff > 1.0

    ;;Electron temperature
    plasma.te = plasma.te > 0. ;[keV]

    ;;Ion temperature
    plasma.ti = plasma.ti > 0. ;[keV]

    if abs(plasma.time - inp.time) gt 0.02 then begin
        warn,'Plasma time and input time do not match'
        print,'Input time: ',inp.time
        print,'Plasma time: ',plasma.time
    endif

    plasma = create_struct(plasma, grid)

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid plasma parameters. Exiting...'
    endif else begin
        success,'Plasma parameters are valid'
    endelse
END

