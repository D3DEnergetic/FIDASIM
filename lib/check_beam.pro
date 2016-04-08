PRO check_beam, inp, nbi, err_status
    ;+#check_beam
    ;+Checks if neutral beam geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **nbi**: neutral beam geometry structure
    ;+ 
    ;+##Output Arguments
    ;+     **err**: error code
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_beam, inputs, nbi, err
    ;+```
    err_status = 0
    info,'Checking beam geometry...'

    zero_string = {dims:0, type:'STRING'}
    zero_int = {dims:0, type:'INT'}
    zero_double = {dims:0, type:'DOUBLE'}
    three_double = {dims:[3], type:'DOUBLE'}
    schema = {data_source:zero_string, $
              name:zero_string, $
              src:three_double,  axis:three_double, $
              divy:three_double, divz:three_double, $
              focy:zero_double,  focz:zero_double, $
              widz:zero_double,  widy:zero_double, $
              shape:zero_int }

    check_struct_schema,schema,nbi,err_status, desc="beam geometry"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    if abs(total(nbi.axis^2.0) - 1.0) gt 1d-5 then begin
        error,'Invalid source axis. Expected norm(axis) == 1'
        err_status = 1
    endif

    if nbi.focz le 0.0 then begin
        error,'focz cannot be in the range (-Inf,0.0]'
        err_status = 1
    endif
   
    if nbi.focy le 0.0 then begin
        error,'focy cannot be in the range (-Inf,0.0]'
        err_status = 1
    endif

    if nbi.shape gt 2 then begin
        error,'Invalid source shape. Expected 1 (rectagular) or 2 (circular)'
        err_status = 1
    endif

    if nbi.widz lt 0 then begin
        error,'Invalid widz. Expected widz > 0'
        err_status = 1
    endif

    if nbi.widy lt 0 then begin
        error,'Invalid widy. Expected widy > 0'
        err_status = 1
    endif

    origin = inp.origin
    uvw_src = nbi.src
    uvw_axis = nbi.axis
    uvw_pos = uvw_src + 200*uvw_axis

    xyz_src = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_src, origin=origin)
    xyz_axis = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_axis)
    xyz_pos = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_pos, origin=origin)
    xyz_center = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,origin, origin=origin)

    dis = sqrt(total((xyz_src - xyz_pos)^2.0))
    BETA=double(asin((xyz_src[2]-xyz_pos[2])/dis))
    ALPHA=double(atan((xyz_pos[1]-xyz_src[1]),(xyz_pos[0]-xyz_src[0])))

    print,'Beam injection start point in machine coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_src
    print,'Point 2m along beam centerline in machine coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_pos
    print,'Machine center in beam grid coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_center
    print,'Beam injection start point in beam grid coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_src
    print,'Point 2m along beam centerline in beam grid coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_pos
    print,'Beam grid rotation angles that would align it with the beam centerline'
    print, ALPHA/!DPI*180,FORMAT='("    alpha = ",F14.10,"°")'
    print, BETA/!DPI*180,FORMAT='("    beta = ",F14.10,"°")'

    ;; Calculate grid center rc and sides length dr
    dr = [inp.xmax-inp.xmin,inp.ymax-inp.ymin,inp.zmax-inp.zmin]
    rc = [inp.xmin,inp.ymin,inp.zmin] + 0.5*dr

    ;; Check if beam centerline intersects beam grid
    aabb_intersect,rc,dr,xyz_src,xyz_axis,length,r_enter,r_exit

    print,'Beam centerline - grid intersection length'
    print,f='("    length = ",F8.3)',length
    if length le 10.0 then begin
        error,'Beam centerline does not intersect grid'
        err_status = 1
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid beam geometry. Exiting...'
    endif else begin
        success,'Beam geometry is valid'
    endelse

END

