PRO check_neutron_collimator, inp, nc
    ;+#check_neutron_collimator
    ;+Checks if Neutron Collimator geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **nc**: Neutron Collimator geometry structure
    ;+ 
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_neutron_collimator, inputs, nc
    ;+```

    err_status=0
    info,'Checking Neutron Collimator geometry...'

    w = where("nchan" eq strlowcase(TAG_NAMES(nc)),nw)
    if nw eq 0 then begin
        error,'"nchan" is missing from the Neutron Collimator geometry'
        err_status = 1
        goto, GET_OUT
    endif
    
    nchan = nc.nchan
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

    check_struct_schema,schema,nc,err_status, desc="Neutron Collimator geometry"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    ;; Check detector/aperture shape
    w = where(nc.d_shape gt 2 or nc.d_shape eq 0,nw)
    if nw ne 0 then begin
        error,'Invalid detector shape. Expected 1 (rectagular) or 2 (circular)'
        print,'Invalid indices: ', w
        err_status = 1
    endif

    w = where(nc.a_shape gt 2 or nc.a_shape eq 0,nw)
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
        uvw_det = nc.d_cent[*,i]
        d_e1 = nc.d_redge[*,i] - uvw_det
        d_e2 = nc.d_tedge[*,i] - uvw_det

        uvw_aper = nc.a_cent[*,i]
        a_e1 = nc.a_redge[*,i] - uvw_aper
        a_e2 = nc.a_tedge[*,i] - uvw_aper

        uvw_dir = uvw_aper - uvw_det

        ;;Rotate chords into beam grid coordinates
        xyz_aper = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_aper,origin=inp.origin)
        xyz_det  = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_det,origin=inp.origin)
        xyz_dir = xyz_aper - xyz_det
        xyz_dir = xyz_dir/sqrt(total(xyz_dir*xyz_dir))

        ;; Check if nc chord intersects beam grid
        aabb_intersect,rc,dr,xyz_det,xyz_dir,length,r_enter,r_exit
        if length le 0.0 then begin
            err_arr[i] = 1
        endif

        ;; Check if Neutron Collimator detector is pointing in the right direction
        d_enter = sqrt(total((r_enter - xyz_aper)^2))
        d_exit = sqrt(total((r_exit - xyz_aper)^2))
        if d_exit lt d_enter then begin
            err_arr[i] = 1
        endif

        ;; Check that the detector and aperture point in the right direction
        d_e3 = crossp(d_e1,d_e2)
        a_e3 = crossp(a_e1,a_e2)
        d_dp = total(uvw_dir*d_e3)
        a_dp = total(uvw_dir*a_e3)
        dp = total(d_e3*a_e3)
        if (a_dp le 0.0) or (d_dp le 0.0) or (dp le 0.0) then begin
            error,'The detector and/or aperture plane normal vectors are pointing in the wrong direction. The Neutron Collimator definition is incorrect.'
            err_arr[i] = 1
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
        print,'    ',nc.id[ww]
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid Neutron Collimator geometry. Exiting...',/halt
    endif else begin
        success,'Neutron Collimator geometry is valid'
    endelse
END

