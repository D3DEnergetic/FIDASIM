PRO check_distribution, inp, grid, dist
    ;+#check_distribution
    ;+Checks if distribution structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: Input structure
    ;+ 
    ;+     **grid**: Interpolation grid structure
    ;+ 
    ;+     **dist**: Fast-ion distribution structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_distribution, inputs, grid, dist
    ;+```

    err_status = 0
    info,'Checking fast-ion distribution'

    w = where("type" eq strlowcase(TAG_names(dist)),nw)
    if nw eq 0 then begin
        error,'"type" is missing from the fast-ion distribution'
        err_status = 1
        goto, GET_OUT
    endif
    dist_type = dist.type

    CASE dist_type OF
        1: BEGIN
            print, 'Using a Guiding Center Fast-ion Density Function'
            w = where("nenergy" eq strlowcase(TAG_names(dist)),nw)
            if nw eq 0 then begin
                error,'"nenergy" is missing from the fast-ion distribution'
                err_status = 1
                goto, GET_OUT
            endif

            w = where("npitch" eq strlowcase(TAG_names(dist)),nw)
            if nw eq 0 then begin
                error,'"npitch" is missing from the fast-ion distribution'
                err_status = 1
                goto, GET_OUT
            endif
            np =  dist.npitch
            nen = dist.nenergy
            nr = grid.nr
            nz = grid.nz
            nphi = grid.nphi
            zero_string = {dims:0, type:'STRING'}
            zero_int = {dims:0, type:'INT'}
            zero_double = {dims:0, type:'DOUBLE'}
            nrnznphi_double = {dims:[nr,nz,nphi], type:'DOUBLE'}
            schema = {type:zero_int, $
                      nenergy:zero_int, $
                      npitch:zero_int, $
                      energy:{dims:[nen], type:'DOUBLE'},$
                      pitch:{dims:[np], type:'DOUBLE'}, $
                      denf:nrnznphi_double, $ 
                      f:{dims:[nen,np,nr,nz,nphi], type:'DOUBLE'}, $
                      time:zero_double, $
                      data_source:zero_string}

            check_struct_schema,schema,dist,err_status, desc="fast-ion distribution"
            if err_status eq 1 then begin
                goto, GET_OUT
            endif
            dist = create_struct(dist, grid)
        END
        2: BEGIN
            print, 'Using Guiding Center Monte Carlo fast-ion distribution'
            w = where("nparticle" eq strlowcase(TAG_names(dist)),nw)
            if nw eq 0 then begin
                error,'"nparticle" is missing from the fast-ion distribution'
                err_status = 1
                goto, GET_OUT
            endif

            npart = dist.nparticle
            zero_int = {dims:0, type:'INT'}
            zero_long = {dims:0, type:'LONG'}
            zero_string = {dims:0, type:'STRING'}
            zero_double = {dims:0, type:'DOUBLE'}
            npart_double = {dims:[npart], type:'DOUBLE'}
            npart_int = {dims:[npart], type:'INT'}
            schema = {type:zero_int, $
                      nparticle:zero_long, $
                      nclass:zero_int, $
                      time:zero_double, $
                      energy:npart_double, $
                      pitch:npart_double, $
                      r:npart_double, $
                      z:npart_double, $
                      weight:npart_double, $
                      class:npart_int,$ 
                      data_source:zero_string}

            check_struct_schema,schema,dist,err_status, desc="fast-ion distribution"
            if err_status eq 1 then begin
                goto, GET_OUT
            endif
            print,'Number of MC particles: ',npart
        END
        3: BEGIN
            print, 'Using Full Orbit Monte Carlo fast-ion distribution'
            w = where("nparticle" eq strlowcase(TAG_names(dist)),nw)
            if nw eq 0 then begin
                error,'"nparticle" is missing from the fast-ion distribution'
                err_status = 1
                goto, GET_OUT
            endif

            npart = dist.nparticle
            zero_int = {dims:0, type:'INT'}
            zero_long = {dims:0, type:'LONG'}
            zero_string = {dims:0, type:'STRING'}
            zero_double = {dims:0, type:'DOUBLE'}
            npart_double = {dims:[npart], type:'DOUBLE'}
            npart_int = {dims:[npart], type:'INT'}
            schema = {type:zero_int, $
                      nparticle:zero_long, $
                      nclass:zero_int, $
                      time:zero_double, $
                      vr:npart_double, $
                      vt:npart_double, $
                      vz:npart_double, $
                      r:npart_double, $
                      z:npart_double, $
                      weight:npart_double, $
                      class:npart_int, $ 
                      data_source:zero_string}

            check_struct_schema,schema,dist,err_status, desc="fast-ion distribution"
            if err_status eq 1 then begin
                goto, GET_OUT
            endif
            print,'Number of MC particles: ',npart
        END
        ELSE: BEGIN
            error,'Invalid distribution type. Expected '+ $
                  '1 (Guiding Center Density Function), '+ $
                  '2 (Guiding Center Monte Carlo), or '+ $
                  '3 (Full Orbit Monte Carlo)'
            err_status = 1
            goto, GET_OUT
            END
        ENDCASE

    if dist.data_source eq '' then begin
        error, 'Invalid data source. An empty string is not a data source.'
        err_status = 1
    endif
 
    if abs(dist.time - inp.time) gt 0.02 then begin
        warn,'Distribution time and input time do not match'
        print,'Input time: ',inp.time
        print,'Distribution time: ',dist.time
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid fast-ion distribution. Exiting...',/halt
    endif else begin
        success,'Fast-ion distribution is valid'
    endelse
END

