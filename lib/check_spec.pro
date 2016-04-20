PRO check_spec, inp, chords, err_status
    ;+#check_spec
    ;+Check if spectral geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **chords**: spectral geometry structure
    ;+ 
    ;+##Output Arguments
    ;+     **err**: error code
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_spec, inputs, chords, err
    ;+```
    err_status = 0
    info,'Checking FIDA/BES inputs...'

    w = where("nchan" eq strlowcase(TAG_NAMES(chords)),nw)
    if nw eq 0 then begin
        error,'"nchan" is missing from the FIDA/BES geometry'
        err_status = 1
        goto, GET_OUT
    endif

    w = where("system" eq strlowcase(TAG_NAMES(chords)),nw)
    if nw eq 0 then begin
        error,'"system" is missing from the FIDA/BES geometry'
        err_status = 1
        goto, GET_OUT
    endif

    nchan = chords.nchan
    nsys = size(chords.system,/dim)
    zero_string = {dims:0, type:'STRING'}
    zero_long = {dims:0, type:'LONG'}
    nchan_double = {dims:[nchan], type:'DOUBLE'}
    schema = {data_source:zero_string, $
              nchan:zero_long, $
              system:{dims:nsys,type:'STRING'}, $
              lens:{dims:[3,nchan], type:'DOUBLE'}, $
              axis:{dims:[3,nchan], type:'DOUBLE'}, $
              sigma_pi:nchan_double, $
              spot_size:nchan_double, $
              radius:nchan_double}

    check_struct_schema,schema,chords,err_status, desc="FIDA/BES geometry"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    err_arr=intarr(nchan)
    cross_arr=intarr(nchan)
    uvw_lens = chords.lens
    uvw_axis = chords.axis

    ;;ROTATE CHORDS INTO BEAM GRID COORDINATES
    xyz_lens = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_lens,origin=inp.origin)
    xyz_axis = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_axis)

    ;; Calculate grid center rc and sides length dr
    dr = [inp.xmax-inp.xmin,inp.ymax-inp.ymin,inp.zmax-inp.zmin]
    rc = [inp.xmin,inp.ymin,inp.zmin] + 0.5*dr

    for i=0,nchan-1 do begin
        chan_str = strcompress(string(i),/remove_all)
        if abs(total(uvw_axis[*,i]^2.0) - 1.0) gt 1d-5 then begin
            error,'Invalid optical axis at index '+chan_str+'. Expected norm(axis) == 1'
            print, total(uvw_axis[*,i]^2.0) - 1.0
            err_arr[i] = 1
        endif

        ;; Check if viewing chord intersects beam grid
        aabb_intersect,rc,dr,xyz_lens[*,i],xyz_axis[*,i],length,r_enter,r_exit
        if length le 0.0 then begin
            warn,'Chord at index '+chan_str+' does not cross the beam grid'
            cross_arr[i] = 1
        endif
    endfor

    w = where(cross_arr eq 0.0,nw,complement=ww,ncomplement=nww)
    print,f='(i3," out of ",i3," chords crossed the beam grid")',nw,nchan
    if nw eq 0 then begin
        error,'No channels intersect the beam grid'
        err_status = 1
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid FIDA/BES geometry. Exiting...'
    endif else begin
        success,'FIDA/BES geometry is valid'
    endelse
END

