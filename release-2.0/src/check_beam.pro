PRO check_beam, inp, nbi
    ;+#check_beam
    ;+Checks if neutral beam geometry structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **nbi**: neutral beam geometry structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_beam, inputs, nbi
    ;+```
    err_status = 0
    info,'Checking beam geometry...'

    na = nbi.naperture
    zero_string = {dims:0, type:'STRING'}
    zero_int = {dims:0, type:'INT'}
    zero_double = {dims:0, type:'DOUBLE'}
    three_double = {dims:[3], type:'DOUBLE'}
    na_double = {dims:[na], type:'DOUBLE'}
    na_int = {dims:[na], type:'INT'}
    schema = {data_source:zero_string, $
              name:zero_string, shape:zero_int, $
              src:three_double,  axis:three_double, $
              divy:three_double, divz:three_double, $
              focy:zero_double,  focz:zero_double, $
              widz:zero_double,  widy:zero_double, $
              naperture:zero_int, ashape:na_int, $
              awidy:na_double, awidz:na_double, $
              aoffy:na_double, aoffz:na_double, $
              adist:na_double }

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

    if nbi.shape gt 2 or nbi.shape eq 0 then begin
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

    w = where(nbi.ashape gt 2 or nbi.ashape eq 0,nw)
    if nw gt 0 then begin
        error,'Invalid aperture shape. Expected 1 (rectangular) or 2 (circular)'
        err_status = 1
    endif

    w = where(nbi.awidy lt 0, nw)
    if nw gt 0 then begin
        error, 'Invalid awidy. Expected awidy >= 0.0'
        err_status = 1
    endif

    w = where(nbi.awidz lt 0, nw)
    if nw gt 0 then begin
        error, 'Invalid awidz. Expected awidz >= 0.0'
        err_status = 1
    endif

    origin = inp.origin
    uvw_src = nbi.src
    uvw_axis = nbi.axis
    if nbi.naperture eq 0 then begin
        uvw_pos = uvw_src + 100*uvw_axis
    endif else begin
        uvw_pos = uvw_src + nbi.adist[0]*uvw_axis
    endelse

    xyz_src = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_src, origin=origin)
    xyz_axis = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_axis)
    xyz_pos = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,uvw_pos, origin=origin)
    xyz_center = uvw_to_xyz(inp.alpha,inp.beta,inp.gamma,[0.0,0.0,0.0], origin=origin)

    dis = sqrt(total((xyz_src - xyz_pos)^2.0))
    BETA=double(asin((xyz_src[2]-xyz_pos[2])/dis))
    ALPHA=double(atan((xyz_pos[1]-xyz_src[1]),(xyz_pos[0]-xyz_src[0])))

    print,'Machine center in beam grid coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_center
    print,'Beam injection start point in machine coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_src
    print,'Beam injection start point in beam grid coordinates'
    print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_src

    if nbi.naperture ne 0 then begin
        print,'First aperture position in machine coordinates'
        print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_pos
        print,'First aperture position in beam grid coordinates'
        print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_pos
    endif else begin
        print,'Position of point 100cm along beam centerline in machine coordinates'
        print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', uvw_pos
        print,'Position of point 100cm along beam centerline in beam grid coordinates'
        print, f='("    [",F9.3,",",F9.3,",",F9.3,"]")', xyz_pos
    endelse

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
        error,'Invalid beam geometry. Exiting...',/halt
    endif else begin
        success,'Beam geometry is valid'
    endelse

END

