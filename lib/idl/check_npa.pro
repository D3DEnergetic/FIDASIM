PRO check_npa, inp, npa
    ;+#check_npa
    ;+Checks if NPA geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **npa**: NPA geometry structure
    ;+ 
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_npa, inputs, npa
    ;+```

    err_status=0
    info,'Checking NPA geometry...'

    w = where("nchan" eq strlowcase(TAG_NAMES(npa)),nw)
    if nw eq 0 then begin
        error,'"nchan" is missing from the NPA geometry'
        err_status = 1
        goto, GET_OUT
    endif
    
    nchan = npa.nchan
    zero_string = {dims:0, type:'STRING'}
    zero_long = {dims:0, type:'LONG'}
    schema = {data_source:zero_string, $
              nchan:zero_long, $
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
              radius:{dims:[nchan], type:'DOUBLE'} }

    check_struct_schema,schema,npa,err_status, desc="NPA geometry"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    ;; Check detector/aperture shape
    w = where(npa.d_shape gt 2 or npa.d_shape eq 0,nw)
    if nw ne 0 then begin
        error,'Invalid detector shape. Expected 1 (rectagular) or 2 (circular)'
        print,'Invalid indices: ', w
        err_status = 1
    endif

    w = where(npa.a_shape gt 2 or npa.a_shape eq 0,nw)
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
        uvw_det = npa.d_cent[*,i]
        d_e1 = npa.d_redge[*,i] - uvw_det
        d_e2 = npa.d_tedge[*,i] - uvw_det

        uvw_aper = npa.a_cent[*,i]
        a_e1 = npa.a_redge[*,i] - uvw_aper
        a_e2 = npa.a_tedge[*,i] - uvw_aper

        ;;Rotate chords into beam grid coordinates
        xyz_aper = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_aper,origin=inp.origin)
        xyz_det  = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_det,origin=inp.origin)
        xyz_dir = xyz_aper - xyz_det
        xyz_dir = xyz_dir/sqrt(total(xyz_dir*xyz_dir))

        ;; Check if npa chord intersects beam grid
        aabb_intersect,rc,dr,xyz_det,xyz_dir,length,r_enter,r_exit
        if length le 0.0 then begin
            err_arr[i] = 1
        endif

        ;; Check that the detector and aperture point in the same direction
        d_e3 = crossp(d_e1,d_e2)
        a_e3 = crossp(a_e1,a_e2)
        dp = total(d_e3*a_e3)
        if dp le 0.0 then begin
           warn,'The dot product of the detector and aperture plane normal vectors is negative. The NPA definition may be incorrect.'
        endif
    endfor

    w = where(err_arr eq 0.0,nw,complement=ww,ncomplement=nww)
    print,f='(i3," out of ",i3," channels crossed the beam grid")',nw,nchan
    if nw eq 0 then begin
        error,'No channels intersect the beam grid'
        err_status = 1
    endif

    if nww gt 0 then begin
        warn,'Some channels did not intersect the beam grid'
        print,'Number missed: ',nww
        print,'Missed channels:'
        print,'    ',npa.id[ww]
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid NPA geometry. Exiting...',/halt
    endif else begin
        success,'NPA geometry is valid'
    endelse
END

