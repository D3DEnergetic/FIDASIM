PRO check_cfpd, inp, cfpd
    ;+#check_cfpd
    ;+Checks if cfpd geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **cfpd**: cfpd geometry structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_cfpd, inputs, cfpd

    err_status=0
    info,'Checking cfpd geometry...'

    w = where("nchan" eq strlowcase(TAG_NAMES(cfpd)),nw)
    if nw eq 0 then begin
        error,'"nchan" is missing from the cfpd geometry'
        err_status = 1
        goto, GET_OUT
    endif

    nchan = cfpd.nchan
    nrays = cfpd.nrays
    nsteps = cfpd.nsteps
    nenergy = cfpd.nenergy
    zero_string = {dims:0, type:'STRING'}
    zero_long = {dims:0, type:'LONG'}
    zero_double = {dims:0,type:'DOUBLE'}
    schema = {data_source:zero_string, $
              nchan:zero_long, $
              nrays:zero_long, $
              nsteps:zero_long, $
              nenergy:zero_long, $
              dl:zero_double, $
              system:zero_string, $
              id:{dims:[nchan], type:'STRING'}, $
              a_shape:{dims:[nchan], type:'INT'},$
              d_shape:{dims:[nchan], type:'INT'}, $
              a_tedge:{dims:[3,nchan], type:'DOUBLE'}, $
              a_redge:{dims:[3,nchan], type:'DOUBLE'}, $
              a_cent:{dims:[3,nchan], type:'DOUBLE'}, $
              d_tedge:{dims:[3,nchan], type:'DOUBLE'}, $
              d_redge:{dims:[3,nchan], type:'DOUBLE'}, $
              d_cent:{dims:[3,nchan], type:'DOUBLE'}, $
              radius:{dims:[nchan], type:'DOUBLE'}, $
              sightline:{dims:[nenergy,6,nsteps,nrays,nchan], type:'DOUBLE'}, $
              daomega:{dims:[nenergy,nrays,nchan], type:'DOUBLE'}, $
              nactual:{dims:[nenergy,nrays,nchan], type:'LONG'}, $
              earray:{dims:[nenergy], type:'DOUBLE'} $
              }

    check_struct_schema,schema,cfpd,err_status, desc="cfpd geometry"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    ;; Check detector/aperture shape
    w = where(cfpd.d_shape gt 2 or cfpd.d_shape eq 0,nw)
    if nw ne 0 then begin
        error,'Invalid detector shape. Expected 1 (rectagular) or 2 (circular)'
        print,'Invalid indices: ', w
        err_status = 1
    endif

    w = where(cfpd.a_shape gt 2 or cfpd.a_shape eq 0,nw)
    if nw ne 0 then begin
        error,'Invalid aperture shape. Expected 1 (rectagular) or 2 (circular)'
        print,'Invalid indices: ', w
        err_status = 1
    endif

    ;; Calculate grid center rc and sides length dr
    dr = [inp.xmax-inp.xmin,inp.ymax-inp.ymin,inp.zmax-inp.zmin]
    rc = [inp.xmin,inp.ymin,inp.zmin] + 0.5*dr
    err_arr = dblarr(nchan)
    for i=0, nchan-1 do begin
        uvw_det = cfpd.d_cent[*,i]
        d_e1 = cfpd.d_redge[*,i] - uvw_det
        d_e2 = cfpd.d_tedge[*,i] - uvw_det

        uvw_aper = cfpd.a_cent[*,i]
        a_e1 = cfpd.a_redge[*,i] - uvw_aper
        a_e2 = cfpd.a_tedge[*,i] - uvw_aper

        uvw_dir = uvw_aper - uvw_det

        ;; Check that the detector and aperture point in the right direction
        d_e3 = crossp(d_e1,d_e2)
        a_e3 = crossp(a_e1,a_e2)
        d_dp = total(uvw_dir*d_e3)
        a_dp = total(uvw_dir*a_e3)
        dp = total(d_e3*a_e3)
        if (a_dp le 0.0) or (d_dp le 0.0) or (dp le 0.0) then begin
            error,'The detector and/or aperture plane normal vectors are pointing in the wrong direction. The cfpd definition is incorrect.'
            err_arr[i] = 1
        endif
    endfor

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid cfpd geometry. Exiting...',/halt
    endif else begin
        success,'cfpd geometry is valid'
    endelse
END

