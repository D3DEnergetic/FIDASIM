;+ This file contains PREFIDA, a procedure that checks and creates FIDASIM {!../VERSION!} inputs
;+
FUNCTION colored, str, c=c, s=s
    ;+***
    ;+##`colored(str, c="w", s="n")`
    ;+Creates colored string
    ;+###Arguments
    ;+     **str**: String to be colored
    ;+
    ;+###Keyword Arguments
    ;+     **c**: Foreground color code
    ;+
    ;+     **s**: Style code
    ;+
    ;+####Foreground Color Codes
    ;+     **k**: Black,
    ;+     **r**: Red,
    ;+     **g**: Green,
    ;+     **y**: Yellow,
    ;+     **b**: Blue,
    ;+     **m**: Magenta,
    ;+     **c**: Cyan,
    ;+     **w**: White
    ;+####Style Format Codes
    ;+     **n**: Normal, 
    ;+     **b**: Bright, 
    ;+     **d**: Dim, 
    ;+     **i**: Italics, 
    ;+     **u**: Underline, 
    ;+     **r**: Reverse, 
    ;+     **h**: Hidden, 
    ;+     **s**: Strikethrough 
    ;+###Return Value
    ;+     Colored string
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> blue_bright_hello = colored("Hello",c="b",s="b")
    ;+```
    if not keyword_set(c) then c='w' ;;Foreground Color
    if not keyword_set(s) then s='n' ;;Style
    esc=string(27b)
    back=esc+"[0m"

    style = {n:'0',b:'1',d:'2',i:'3',u:'4',r:'7',h:'8',s:'9'}
    sTags = ["n","b","d","i","u","r","h","s"]

    fgColors = { k:'30',r:'31',g:'32',y:'33',$
                b:'34',m:'35',c:'36',w:'37'}
    fgTags = ["k","r","g","y","b","m","c","w"]

    sIndex = where(s eq sTags)
    fgIndex = where(c eq fgTags)
    if sIndex eq -1 then sIndex=0
    if fgIndex eq -1 then fgIndex=7

    return, esc+"["+style.(sIndex)+";"+fgColors.(fgIndex)+"m"+str+back
END

PRO success, str
    ;+***
    ;+##`success, str`
    ;+Print a success message
    ;+###Arguments
    ;+     **str**: message
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> success, "Yay!!!"
    ;+```
    print, colored('SUCCESS: '+str,c='g')
END

PRO warn, str
    ;+***
    ;+##`warn, str`
    ;+Print a warning message
    ;+###Arguments
    ;+     **str**: message
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> warn, "This may be a problem"
    ;+```
    print, colored('WARNING: '+str,c='y')
END

PRO error, str
    ;+***
    ;+##`error, str`
    ;+Print a error message
    ;+###Arguments
    ;+     **str**: message
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> error, "=("
    ;+```
    print, colored('ERROR: '+str,c='r')
END

PRO info, str
    ;+***
    ;+##`info, str`
    ;+Print a informational message
    ;+###Arguments
    ;+     **str**: message
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> info, "This is an informative message"
    ;+```
    print, colored('INFO: '+str,c='b')
END

PRO check_struct_schema, schema, s, err_status, desc=desc
    ;+***
    ;+##`check_struct_schema, schema, s, err, desc="structure"`
    ;+ Check structure `s` is formatted according to `schema`
    ;+###Input Arguments
    ;+     **schema**: structure schema
    ;+
    ;+     **s**: structure to check
    ;+
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Keyword Arguments
    ;+     **desc**: description of structure `s`
    ;+
    ;+ Example usage:
    ;+```idl
    ;+IDL> s = {a:0, b:[1.d0,2.d0], c:"example"}
    ;+IDL> schema = {a:{dims:0,type:"INT"}, b:{dims:[2],type:"DOUBLE"}, c:{dims:0,type:"STRING"}  }
    ;+
    ;+IDL> check_struct_schema, schema, s, err, desc="Example structure"
    ;+IDL> print, err
    ;+    0
    ;+```
    if not keyword_set(struct_name) then desc = 'structure'

    err_status = 0
    schema_tags = strlowcase(TAG_NAMES(schema)) 
    tags = strlowcase(TAG_NAMES(s)) 

    for i=0,n_elements(tags)-1 do begin
        w=where(tags[i] eq schema_tags,nw)
        if nw eq 0 then begin
            info,'Extra variable "'+tags[i]+'" found in '+desc
        endif
    endfor

    for i=0,n_elements(schema_tags)-1 do begin
        w=where(schema_tags[i] eq tags,nw)
        if nw eq 0 then begin
            error,'"'+schema_tags[i]+'" is missing from the '+desc
            err_status = 1
        endif else begin
            ;; Check dimensions
            ww = where((size(s.(w),/dim) eq schema.(i).dims) ne 1,nww) 
            if nww ne 0 then begin
                error,'"'+schema_tags[i]+'" has the wrong dimensions. Expected ('+ $
                      strjoin(strcompress(string(schema.(i).dims),/remove_all),',')+')'
                print,'size('+schema_tags[i]+') = ',size(s.(w),/dim)
                err_status = 1
            endif
            ;; Check type
            tname = size(s.(w),/tname)
            if tname ne schema.(i).type then begin
                error,'"'+schema_tags[i]+'" has the wrong type. Expected '+schema.(i).type
                print,'type('+schema_tags[i]+') = '+tname
                err_status = 1
            endif
            ;; Check for NaNs or Inf
            if tname ne 'STRING' and tname ne 'STRUCT' then begin
                ww = where(finite(s.(w)) eq 0,nww) 
            endif else nww = 0
            if nww ne 0 then begin
                error,'NaN or Infinity detected in "'+schema_tags[i]+'"'
                err_status = 1
            endif
        endelse
    endfor

END

FUNCTION tb_zyx, a, b, g
    ;+***
    ;+##`tb_zyx(alpha, beta, gamma)`
    ;+Calculates Tait-Bryan z-y'-x" active rotation matrix given rotation angles `alpha`,`beta`,`gamma` in radians
    ;+###Arguments
    ;+     **alpha**: rotation angle about z [radians]
    ;+
    ;+     **beta**: rotation angle about y' [radians]
    ;+
    ;+     **gamma**: rotation angle about x" [radians]
    ;+
    ;+###Return value
    ;+     Rotation Matrix
    ;+
    ;+###Example Usage
    ;+```idl
    ;+ IDL> rot_mat = tb_zyx(!DPI/2, 0.0, !DPI/3)
    ;+```
    sa = sin(a) & ca = cos(a)
    sb = sin(b) & cb = cos(b)
    sg = sin(g) & cg = cos(g)

    R = dblarr(3,3)
    R[0,0] = ca*cb & R[1,0] = ca*sb*sg - cg*sa & R[2,0] = sa*sg + ca*cg*sb
    R[0,1] = cb*sa & R[1,1] = ca*cg + sa*sb*sg & R[2,1] = cg*sa*sb - ca*sg
    R[0,2] = -sb   & R[1,2] = cb*sg            & R[2,2] = cb*cg

    return, R
END

FUNCTION tile_array, arr, ncol, nrow
    ;+***
    ;+##`tile_array(arr, ncol, nrow)`
    ;+Creates a tiled matrix out of an array or matrix
    ;+###Arguments
    ;+    **arr**: Array or Matrix of size (nx,ny) to be tiled
    ;+
    ;+    **ncol**: Number of columns in the tile
    ;+
    ;+    **nrow**: Number of rows in the tile
    ;+
    ;+###Return Value
    ;+    Tiled array of size (ncol*nx,nrow*ny)
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> arr = indgen(2,2)
    ;+IDL> print, arr
    ;+    0  1
    ;+    2  3
    ;+IDL> tiled_arr = tile_array(arr, 2, 2)
    ;+IDL> print, tiled_arr
    ;+    0  1  0  1
    ;+    2  3  2  3
    ;+    0  1  0  1
    ;+    2  3  2  3
    ;+```
    s = size(arr,/dim)
    if n_elements(s) eq 1 then s = [s,1]
    new_arr = make_array(s[0]*ncol,s[1]*nrow, type=size(arr,/type))
    new_arr[0,0] = arr

    if nrow gt 1 then begin
        for i=1,nrow-1 do begin
            new_arr[0,i*s[1]] = arr
        endfor
    endif

    if ncol gt 1 then begin
        rarr = new_arr[0:s[0]-1,*]
        for j=1,ncol-1 do begin
            new_arr[(j*s[0]),0] = rarr
        endfor
    endif

    return, new_arr

END

FUNCTION xyz_to_uvw, alpha, beta, gamma, xyz, origin = origin
    ;+***
    ;+##`xyz_to_uvw(alpha, beta, gamma, xyz, origin=[0,0,0])`
    ;+Express rotated coordinate `xyz` in non-rotated `uvw` coordinates
    ;+###Arguments
    ;+     **alpha**: Rotation angle about z [radians]
    ;+
    ;+     **beta**: Rotation angle about y' [radians]
    ;+
    ;+     **gamma**: Rotation angle about x" [radians]
    ;+
    ;+     **xyz**: Point in rotated coordinate system
    ;+
    ;+###Keyword Arguments
    ;+     **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates.
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> uvw = xyz_to_uvw(!DPI/2,0.0,!DPI/3,xyz)
    ;+```
    if not keyword_set(origin) then origin=[0.0,0.0,0.0]
    s = size(xyz,/dim)
    if n_elements(s) ne 2 then s=[s,1]
    xyz = transpose(xyz) ;Column Vector

    R = tb_zyx(alpha,beta,gamma)

    uvw = R##xyz 

    return, tranpose(uvw) + tile_array(origin,1,s[1])
END

FUNCTION uvw_to_xyz, alpha, beta, gamma, uvw, origin=origin
    ;+***
    ;+##`uvw_to_xyz(alpha, beta, gamma, uvw, origin=[0,0,0])`
    ;+ Express non-rotated coordinate `uvw` in rotated `xyz` coordinates
    ;+###Arguments
    ;+     **alpha**: Rotation angle about z [radians]
    ;+
    ;+     **beta**: Rotation angle about y' [radians]
    ;+
    ;+     **gamma**: Rotation angle about x" [radians]
    ;+
    ;+     **xyz**: Point in rotated coordinate system
    ;+
    ;+###Keyword Arguments
    ;+     **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates.
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> xyz = uvw_to_xyz(!DPI/2,0.0,!DPI/3,uvw)
    ;+```
    if not keyword_set(origin) then origin=[0.0,0.0,0.0]
    s = size(uvw,/dim)
    if n_elements(s) ne 2 then s=[s,1]
    uvw_shifted = transpose(uvw - tile_array(origin,1,s[1])) 

    R = transpose(tb_zyx(alpha,beta,gamma))

    xyz = R##uvw_shifted

    return, transpose(xyz)
END

PRO aabb_intersect, rc, dr, r0, d0, intersect, r_enter, r_exit
    ;+***
    ;+##`aabb_intersect, rc, dr, r0, d0, length, ri, rf`
    ;+Calculates intersection length of a ray and an axis aligned bounding box (AABB)
    ;+###Input Arguments
    ;+     **rc**: Center of AABB
    ;+
    ;+     **dr**: [length, width, height] of AABB
    ;+
    ;+     **r0**: starting point of ray
    ;+
    ;+     **d0**: direction of ray
    ;+
    ;+###Output Arguments
    ;+     **intersect**: Intersection length of ray and AABB
    ;+
    ;+     **ri**: Optional, ray enterence point
    ;+
    ;+     **rf**: Optional, ray exit point
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> aabb_intersect, [0,0,0], [1,1,1], [-1,0,0], [1,0,0], intersect, ri, rf
    ;+IDL> print, intersect
    ;+    1.0
    ;+IDL> print, ri
    ;+    -0.5  0.0  0.0
    ;+IDL> print, rf
    ;+     0.5  0.0  0.0
    ;+```

    v0 = d0/sqrt(total(d0*d0))

    ;; There are 6 sides to a cube/grid
    side_inter = dblarr(6)
    ;; Intersection points of ray with planes defined by grid
    ipnts = dblarr(3,6)

    ;; Find whether ray intersects each side
    for i=0L,5 do begin
        j = fix(floor(i/2))
        ind = where([0,1,2] ne j)
        if abs(v0[j]) gt 0 then begin
            ;; Intersection point with plane
            ipnts[*,i] = r0 + v0*( ( (rc[j] + ( (i mod 2)-0.5)*dr[j] ) - r0[j])/v0[j] )
            ;; Check if point on plane is within grid side
            if abs(ipnts[ind[0],i] - rc[ind[0]]) le 0.5*dr[ind[0]] and $
               abs(ipnts[ind[1],i] - rc[ind[1]]) le 0.5*dr[ind[1]] then side_inter[i]=1
        endif
    endfor

    intersect = 0.0
    r_enter = r0
    r_exit = r0
    w = where(side_inter ne 0,nw)
    if nw ge 2 then begin
        ;;Find two unique intersection points
        nunique = 0
        for i=0,nw-2 do begin
            if total(ipnts[*,w[0]] eq ipnts[*,w[i+1]]) ne 3 then begin
                w = [w[0],w[i+1]]
                nunique = 2
                break
            end
        end

        if nunique eq 2 then begin
            vi = ipnts[*,w[1]]-ipnts[*,w[0]]
            vi = vi/sqrt(total(vi*vi))
            dot_prod = total(v0*vi)
            if dot_prod gt 0.0 then begin
                r_enter = ipnts[*,w[0]]
                r_exit = ipnts[*,w[1]]
            endif else begin
                r_enter = ipnts[*,w[1]]
                r_exit = ipnts[*,w[0]]
            endelse
            ;; Calculate intersection length
            intersect = sqrt(total((r_exit-r_enter)^2.0))
        endif
    endif
END

PRO check_inputs, inputs, err_status
    ;+***
    ;+##`check_inputs, inputs, err`
    ;+Checks if input structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: input structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> check_inputs, inputs, err
    ;+```
    info,'Checking simulation settings...'
    err_status = 0
    zero_string = {dims:0,type:'STRING'}
    zero_int = {dims:0,type:'INT'}
    zero_long = {dims:0,type:'LONG'}
    zero_double = {dims:0,type:'DOUBLE'}
    three_double = {dims:[3],type:'DOUBLE'}
    schema = {comment:zero_string, $
              shot:zero_long, time:zero_double, $
              runid:zero_string, device:zero_string, $ 
              install_dir:zero_string, tables_file:zero_string, result_dir:zero_string, $
              nlambda:zero_int, lambdamin:zero_double, lambdamax:zero_double, $
              nx:zero_int, ny:zero_int, nz:zero_int, $
              alpha:zero_double, beta:zero_double, gamma:zero_double, $
              origin:three_double, xmin:zero_double, xmax:zero_double, $
              ymin:zero_double, ymax:zero_double, zmin:zero_double, zmax:zero_double, $
              ab:zero_double, ai:zero_double, species_mix:three_double, $
              pinj:zero_double, einj:zero_double, impurity_charge:zero_int, $
              n_fida:zero_long, n_nbi:zero_long, n_dcx:zero_long, $
              n_npa:zero_long, n_halo:zero_long, n_birth:zero_long, $
              ne_wght:zero_int, np_wght:zero_int, nphi_wght:zero_int, $
              emax_wght:zero_double, nlambda_wght:zero_int, $
              lambdamin_wght:zero_double, lambdamax_wght:zero_double, $
              calc_npa:zero_int, calc_fida:zero_int, calc_bes:zero_int, $
              calc_brems:zero_int, calc_birth:zero_int, $
              calc_fida_wght:zero_int, calc_npa_wght:zero_int, $
              dump_dcx:zero_int, load_neutrals:zero_int, verbose:zero_int}

    check_struct_schema, schema, inputs, err_status, desc="simulation settings"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    ;Normalize File Paths
    inputs.result_dir = expand_path(inputs.result_dir)
    inputs.install_dir = expand_path(inputs.install_dir)

    if inputs.alpha gt 2*!DPI or $
       inputs.beta gt 2*!DPI or $
       inputs.gamma gt 2*!DPI then begin
        error,'Angles must be in radians'
        err_status = 1
    endif

    if inputs.lambdamin ge inputs.lambdamax then begin
        error,'Invalid wavelength range. Expected lambdamin < lamdbdamax'
        err_status = 1
    endif

    if inputs.lambdamin_wght ge inputs.lambdamax_wght then begin
        error,'Invalid wavelength range. Expected lambdamin_wght < lamdbdamax_wght'
        err_status = 1
    endif

    if inputs.xmin ge inputs.xmax then begin
        error,'Invalid x range. Expected xmin < xmax'
        err_status = 1
    endif

    if inputs.ymin ge inputs.ymax then begin
        error,'Invalid y range. Expected ymin < ymax'
        err_status = 1
    endif

    if inputs.zmin ge inputs.zmax then begin
        error,'Invalid z range. Expected zmin < zmax'
        err_status = 1
    endif

    if inputs.pinj le 0. or inputs.einj le 0.0 then begin
        error,'The selected source is not on'
        print,'einj = ',inputs.einj
        print,'pinj = ',inputs.pinj
        err_status = 1
    endif

    if abs(total(inputs.species_mix) - 1.0) gt 1.d-3 then begin
        error,'species_mix does not sum to 1.0'
        print,'sum(species_mix) = ',total(inputs.species_mix)
        err_status = 1
    endif

    if inputs.impurity_charge le 1 then begin
        error,'Invalid impurity charge. Expected impurity charge > 1'
        err_status = 1
    endif

    ps = path_sep()
    input_file = inputs.result_dir+ps+inputs.runid+'_inputs.dat'
    equilibrium_file = inputs.result_dir+ps+inputs.runid+'_equilibrium.h5'
    geometry_file = inputs.result_dir+ps+inputs.runid+'_geometry.h5'
    distribution_file = inputs.result_dir+ps+inputs.runid+'_distribution.h5'
    neutrals_file = inputs.result_dir+ps+inputs.runid+'_neutrals.h5'

    inputs = create_struct(inputs,'input_file',input_file, $
                                  'equilibrium_file',equilibrium_file,$ 
                                  'geometry_file',geometry_file, $
                                  'distribution_file',distribution_file, $
                                  'neutrals_file',neutrals_file)

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid simulation settings. Exiting...'
    endif else begin
        success,'Simulation settings are valid'
    endelse

END

PRO check_grid, grid, err_status
    ;+***
    ;+##`check_grid, grid, err`
    ;+Checks if interpolation grid structure is valid
    ;+
    ;+###Input Arguments
    ;+     **grid**: Interpolation grid structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> check_grid, grid, err
    ;+```
    err_status = 0
    info,'Checking interpolation grid...'
  
    w = where("nr" eq strlowcase(TAG_NAMES(grid)),nw)
    if nw eq 0 then begin
        error,'"nr" is missing from the interpolation grid'
        err_status = 1
        goto, GET_OUT
    endif
    
    w = where("nz" eq strlowcase(TAG_NAMES(grid)),nw)
    if nw eq 0 then begin
        error,'"nz" is missing from the interpolation grid'
        err_status = 1
        goto, GET_OUT
    endif

    nr = grid.nr
    nz = grid.nz
    zero_int = {dims:0,type:'INT'}
    schema = {nr:zero_int, nz:zero_int, $
              r2d:{dims:[nr,nz], type:'DOUBLE'}, $
              z2d:{dims:[nr,nz], type:'DOUBLE'}, $
              r:{dims:[nr], type:'DOUBLE'}, $
              z:{dims:[nz], type:'DOUBLE'} }

    check_struct_schema, schema, grid, err_status, desc="interpolation grid"
    if err_status eq 1 then begin
        goto, GET_OUT
    endif

    w = where((indgen(nr) eq sort(grid.r)) ne 1, nw)
    if nw ne 0 then begin
        error,'r is not in ascending order'
        err_status = 1
    endif

    w = where((indgen(nz) eq sort(grid.z)) ne 1, nw)
    if nw ne 0 then begin
        error,'z is not in ascending order'
        err_status = 1
    endif

    w = where((grid.r eq grid.r2d[*,0]) ne 1, nw)
    if nw ne 0 then begin
        error,'r2d is defined incorrectly. Expected r == r2d[*,0]'
        err_status = 1
    endif

    w = where((grid.z eq grid.z2d[0,*]) ne 1, nw)
    if nw ne 0 then begin
        error,'z2d is defined incorrectly. Expected z == z2d[0,*]'
        err_status = 1
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid interpolation grid. Exiting...'
    endif else begin
        success,'Interpolation grid is valid'
    endelse

END

PRO check_beam, inp, nbi, err_status
    ;+***
    ;+##`check_beam, inputs, nbi, err`
    ;+Check if neutral beam geometry structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **nbi**: neutral beam geometry structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
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

PRO check_plasma, inp, grid, plasma, err_status
    ;+***
    ;+##`check_plasma, inputs, grid, plasma, err`
    ;+Checks if plasma paramters structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: Input structure
    ;+
    ;+     **grid**: Interpolation grid structure
    ;+
    ;+     **plasma**: Plasma parameters structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
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

PRO check_fields, inp, grid, fields, err_status
    ;+***
    ;+##`check_fields, inputs, grid, fields, err`
    ;+Check if electromagnetic fields structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: Input structure
    ;+ 
    ;+     **grid**: Interpolation grid structure
    ;+ 
    ;+     **fields**: Electromagnetic fields structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> check_fields, inputs, grid, fields, err
    ;+```

    err_status=0
    info,'Checking electromagnetic fields...'

    nr = grid.nr
    nz = grid.nz
    zero_string = {dims:0, type:'STRING'}
    zero_double = {dims:0, type:'DOUBLE'}
    nrnz_double = {dims:[nr,nz], type:'DOUBLE'}
    nrnz_int = {dims:[nr,nz], type:'INT'}
    schema = {time:zero_double, $
              br:nrnz_double, $
              bt:nrnz_double, $
              bz:nrnz_double, $
              er:nrnz_double, $
              et:nrnz_double, $
              ez:nrnz_double, $
              mask:nrnz_int, $
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
        error,'Invalid electromagnetic fields. Exiting...'
    endif else begin
        success,'Electromagnetic fields are valid'
    endelse
END

PRO check_dist, inp, grid, dist, err_status
    ;+***
    ;+##`check_dist, inputs, grid, dist, err`
    ;+Check if distribution structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: Input structure
    ;+ 
    ;+     **grid**: Interpolation grid structure
    ;+ 
    ;+     **dist**: Fast-ion distribution structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> check_distribution, inputs, grid, dist, err
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
            zero_string = {dims:0, type:'STRING'}
            zero_int = {dims:0, type:'INT'}
            zero_double = {dims:0, type:'DOUBLE'}
            nrnz_double = {dims:[nr,nz], type:'DOUBLE'}
            schema = {type:zero_int, $
                      nenergy:zero_int, $
                      npitch:zero_int, $
                      energy:{dims:[nen], type:'DOUBLE'},$
                      pitch:{dims:[np], type:'DOUBLE'}, $
                      denf:nrnz_double, $ 
                      f:{dims:[nen,np,nr,nz], type:'DOUBLE'}, $
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
        error,'Invalid fast-ion distribution. Exiting...'
    endif else begin
        success,'Fast-ion distribution is valid'
    endelse
END

PRO check_spec, inp, chords, err_status
    ;+***
    ;+##`check_spec, inputs, chords, err`
    ;+Check if spectral geometry structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **chords**: spectral geometry structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
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
    zero_int = {dims:0, type:'INT'}
    nchan_double = {dims:[nchan], type:'DOUBLE'}
    schema = {data_source:zero_string, $
              nchan:zero_int, $
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

PRO check_npa, inp, npa, err_status
    ;+***
    ;+##`check_npa, inputs, npa, err`
    ;+Check if NPA geometry structure is valid
    ;+
    ;+###Input Arguments
    ;+     **inputs**: input structure
    ;+
    ;+     **npa**: NPA geometry structure
    ;+ 
    ;+###Output Arguments
    ;+     **err**: error code
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> check_npa, inputs, npa, err
    ;+```

    err_status=0
    info,'Checking NPA geometry...'

    w = where("nchan" eq strlowcase(TAG_NAMES(npa)),nw)
    if nw eq 0 then begin
        error,'"nchan" is missing from the NPA geometry'
        err_status = 1
        goto, GET_OUT
    endif
    
    w = where("system" eq strlowcase(TAG_NAMES(npa)),nw)
    if nw eq 0 then begin
        error,'"system" is missing from the NPA geometry'
        err_status = 1
        goto, GET_OUT
    endif

    nsys = size(npa.system,/dim)
    nchan = npa.nchan
    zero_string = {dims:0, type:'STRING'}
    zero_int = {dims:0, type:'INT'}
    schema = {data_source:zero_string, $
              nchan:zero_int, $
              system:{dims:nsys, type:'STRING'}, $
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
    w = where(npa.d_shape gt 2,nw)
    if nw ne 0 then begin
        error,'Invalid detector shape. Expected 1 (rectagular) or 2 (circular)'
        print,'Invalid indices: ', w
        err_status = 1
    endif

    w = where(npa.a_shape gt 2,nw)
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
        print,'Missed indices:'
        print,'    ',ww
    endif

    GET_OUT:
    if err_status ne 0 then begin
        error,'Invalid NPA geometry. Exiting...'
    endif else begin
        success,'NPA geometry is valid'
    endelse
END

FUNCTION get_version, fidasim_dir
    ;+***
    ;+##`get_version(fidasim_dir)`
    ;+ Gets FIDASIM version number from git. 
    ;+ Falls back to reading VERSION file when git is not available
    ;+
    ;+###Input Arguments
    ;+    **fidasim_dir**: FIDASIM install directory
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> version = get_version("~/FIDASIM")
    ;+```

    version = ''
    git_dir = fidasim_dir+'/.git'
    spawn,'command -v git ',git_command,/sh
    if file_test(git_command) and file_test(git_dir,/dir) then begin
        spawn,git_command+' --git-dir='+git_dir+' describe --tags --always',version,err_status
    endif else begin
        version_file = fidasim_dir+'/VERSION'
        version = ''
        if file_test(version_file) then begin
            openr, lun, version_file, /get_lun
            readf, lun, version
            free_lun, lun
        endif
    endelse

    return, version

END

PRO write_namelist, filename, inputs
    ;+***
    ;+##`write_namelist, filename, inputs`
    ;+Writes namelist file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the namelist file
    ;+
    ;+     **inputs**: Input structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_namelist, filename, inputs
    ;+```
    info,'Writing namelist file...'

    fidasim_version = get_version(inputs.install_dir)

    openw,55,filename
    printf,55,'!! Created: ', systime()
    printf,55,'!! FIDASIM version: '+fidasim_version
    printf,55,'!! Comment: '+inputs.comment
    printf,55,'&fidasim_inputs'
    printf,55,''
    printf,55,'!! Shot Info'
    printf,55,f='("shot = ", i6 ,"    !! Shot Number")',inputs.shot
    printf,55,f='("time = ", 1f8.5 ,"    !! Time [s]")',inputs.time
    printf,55,"runid = '" + inputs.runid + "'    !! runID"
    printf,55,"result_dir = '" + inputs.result_dir+"'    !! Result Directory"
    printf,55,''
    printf,55,'!! Input Files'
    printf,55,"tables_file = '" + inputs.tables_file+"'    !! Atomic Tables File"
    printf,55,"equilibrium_file = '" + inputs.equilibrium_file +"'    !! File containing plasma parameters and fields"
    printf,55,"geometry_file = '" + inputs.geometry_file +"'    !! File containing NBI and diagnostic geometry"
    printf,55,"distribution_file = '" + inputs.distribution_file +"'    !! File containing fast-ion distribution"
    printf,55,"neutrals_file = '" + inputs.neutrals_file +"'    !! File containing the neutral density"
    printf,55,''
    printf,55,'!! Simulation Switches'
    printf,55,f='("calc_bes = ",i2 , "    !! Calculate Beam Emission Spectra")',inputs.calc_bes
    printf,55,f='("calc_brems = ",i2 , "    !! Calculate Bremsstrahlung")',inputs.calc_brems
    printf,55,f='("calc_fida = ",i2 , "    !! Calculate FIDA Spectra")',inputs.calc_fida
    printf,55,f='("calc_npa = ",i2 , "   !! Calculate NPA")',inputs.calc_npa
    printf,55,f='("calc_birth = ",i2 , "    !! Calculate Birth Profile")',inputs.calc_birth
    printf,55,f='("calc_fida_wght = ",i2 , "    !! Calculate FIDA weights")',inputs.calc_fida_wght
    printf,55,f='("calc_npa_wght = ",i2 , "    !! Calculate NPA weights")',inputs.calc_npa_wght
    printf,55,f='("load_neutrals = ",i2,"    !! Load neutrals from neutrals file")',inputs.load_neutrals
    printf,55,f='("dump_dcx = ",i2,"    !! Dump DCX neutrals and spectra")',inputs.dump_dcx
    printf,55,f='("verbose = ",i2,"    !! Verbose")',inputs.verbose
    printf,55,''
    printf,55,'!! Monte Carlo Settings'
    printf,55,f='("n_fida = ",i9,"    !! Number of FIDA mc particles")',inputs.n_fida
    printf,55,f='("n_npa = ",i9,"    !! Number of NPA mc particles")',inputs.n_npa
    printf,55,f='("n_nbi = ",i9,"    !! Number of NBI mc particles")',inputs.n_nbi
    printf,55,f='("n_halo = ",i9,"    !! Number of HALO mc particles")',inputs.n_halo
    printf,55,f='("n_dcx = ",i9,"     !! Number of DCX mc particles")',inputs.n_dcx
    printf,55,f='("n_birth = ",i9,"    !! Number of BIRTH mc particles")',inputs.n_birth
    printf,55,''
    printf,55,'!! Neutral Beam Settings'
    printf,55,f='("ab = ",1f9.5,"     !! Beam Species mass [amu]")',inputs.ab
    printf,55,f='("pinj = ",1f9.3,"     !! Beam Power [MW]")',inputs.pinj
    printf,55,f='("einj = ",1f9.3,"     !! Beam Energy [keV]")',inputs.einj
    printf,55,f='("species_mix(1) = ",1f9.5,"     !! Beam Species Mix (Full component)")',inputs.species_mix[0]
    printf,55,f='("species_mix(2) = ",1f9.5,"     !! Beam Species Mix (Half component)")',inputs.species_mix[1]
    printf,55,f='("species_mix(3) = ",1f9.5,"     !! Beam Species Mix (Third component)")',inputs.species_mix[2]
    printf,55,''
    printf,55,'!! Plasma Settings'
    printf,55,f='("ai = ",1f9.5,"     !! Ion Species mass [amu]")',inputs.ai
    printf,55,f='("impurity_charge = ",i3,"     !! Impurity Charge")',inputs.impurity_charge
    printf,55,''
    printf,55,'!! Beam Grid Settings'
    printf,55,f='("nx = ",i4,"    !! Number of cells in X direction (Into Plasma)")',inputs.nx
    printf,55,f='("ny = ",i4,"    !! Number of cells in Y direction")',inputs.ny
    printf,55,f='("nz = ",i4,"    !! Number of cells in Z direction")',inputs.nz
    printf,55,f='("xmin = ",1f9.3,"     !! Minimum X value [cm]")',inputs.xmin
    printf,55,f='("xmax = ",1f9.3,"     !! Maximum X value [cm]")',inputs.xmax
    printf,55,f='("ymin = ",1f9.3,"     !! Minimum Y value [cm]")',inputs.ymin
    printf,55,f='("ymax = ",1f9.3,"     !! Maximum Y value [cm]")',inputs.ymax
    printf,55,f='("zmin = ",1f9.3,"     !! Minimum Z value [cm]")',inputs.zmin
    printf,55,f='("zmax = ",1f9.3,"     !! Maximum Z value [cm]")',inputs.zmax
    printf,55,'!! Tait-Bryan Angles for z-y`-x`` rotation'
    printf,55,f='("alpha = ",1f9.5,"     !! Rotation about z-axis [rad]")',inputs.alpha
    printf,55,f='("beta  = ",1f9.5,"     !! Rotation about y`-axis [rad]")',inputs.beta
    printf,55,f='("gamma = ",1f9.5,"     !! Rotation about x``-axis [rad]")',inputs.gamma
    printf,55,'!! Beam Grid origin in machine coordinates (cartesian)'
    printf,55,f='("origin(1) = ",1f9.3,"     !! U value [cm]")',inputs.origin[0]
    printf,55,f='("origin(2) = ",1f9.3,"     !! V value [cm]")',inputs.origin[1]
    printf,55,f='("origin(3) = ",1f9.3,"     !! W value [cm]")',inputs.origin[2]
    printf,55,''
    printf,55,'!! Wavelength Grid Settings'
    printf,55,f='("nlambda = ",1i5,"    !! Number of Wavelengths")',inputs.nlambda
    printf,55,f='("lambdamin = ",1f9.3,"    !! Minimum Wavelength [nm]")',inputs.lambdamin
    printf,55,f='("lambdamax = ",1f9.3,"    !! Maximum Wavelength [nm]")',inputs.lambdamax
    printf,55,''
    printf,55,'!! Weight Function Settings'
    printf,55,f='("ne_wght = ",i9,"    !! Number of Energies for Weights")',inputs.ne_wght
    printf,55,f='("np_wght = ",i9,"    !! Number of Pitches for Weights")',inputs.np_wght
    printf,55,f='("nphi_wght = ",i9,"    !! Number of Gyro-angles for Weights")',inputs.nphi_wght
    printf,55,f='("emax_wght = ",1f9.2,"    !! Maximum Energy for Weights [keV]")',inputs.emax_wght
    printf,55,f='("nlambda_wght = ",1i5,"    !! Number of Wavelengths for Weights ")',$
              inputs.nlambda_wght
    printf,55,f='("lambdamin_wght = ",1f9.3,"    !! Minimum Wavelength for Weights [nm]")',$
              inputs.lambdamin_wght
    printf,55,f='("lambdamax_wght = ",1f9.3,"    !! Maximum Wavelength for Weights [nm]")',$
              inputs.lambdamax_wght
    printf,55,''
    printf,55,'/'
    printf,55,''
    close,55
    success,'Namelist file created: '+filename
END

PRO write_geometry, filename, nbi, spec=spec, npa=npa
    ;+***
    ;+##`write_geometry, filename, nbi, spec=spec, npa=npa`
    ;+Write geometry values to a HDF5 file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the geometry file
    ;+
    ;+     **nbi**: NBI geometry structure
    ;+
    ;+###Keyword Arguments
    ;+     **spec**: Optional, Spectral geometry structure
    ;+
    ;+     **npa**: Optional, NPA geometry structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_geometry, filename, nbi, spec=spec, npa=npa
    ;+```
    info,'Writing geometry file...'

    ;; Create attributes
    root_atts = {attribute,obj:'/', $
                name:'description',$
                data:'Geometric quantities for FIDASIM'}

    ;; NBI attributes
    nbi_desc = {attribute,obj:'/nbi', $
                name:'description',$
                data:'Neutral Beam Geometry'}

    nbi_cs = {attribute,obj:'/nbi', $
              name:'coordinate_system',$
              data:'Right-handed cartesian'}

    nbi_ds_desc = {attribute,obj:'/nbi/data_source', $
                   name:'description', $
                   data:'Source of the NBI geometry'}

    nbi_name_desc = {attribute,obj:'/nbi/name', $
                     name:'description', $
                     data:'Beam name'}

    nbi_src_desc = {attribute,obj:'/nbi/src', $
                    name:'description', $
                    data:'Position of the center of the beam source grid'}
    nbi_src_unit = {attribute,obj:'/nbi/src', $
                    name:'units', $
                    data:'cm'}

    nbi_axis_desc = {attribute,obj:'/nbi/axis', $
                     name:'description', $
                     data:'Axis of the beam centerline: Centerline(t) = src + axis*t '}
    nbi_axis_unit = {attribute,obj:'/nbi/axis', $
                     name:'units', $
                     data:'cm'}

    nbi_focy_desc = {attribute,obj:'/nbi/focy', $
                     name:'description', $
                     data:'Horizonal focal length of the beam'}
    nbi_focy_unit = {attribute,obj:'/nbi/focy', $
                     name:'units', $
                     data:'cm'}

    nbi_focz_desc = {attribute,obj:'/nbi/focz', $
                     name:'description', $
                     data:'Vertical focal length of the beam'}
    nbi_focz_unit = {attribute,obj:'/nbi/focz', $
                     name:'units', $
                     data:'cm'}
     
    nbi_divy_desc = {attribute,obj:'/nbi/divy', $
                     name:'description', $
                     data:'Horizonal divergences of the beam. One for each energy component'}
    nbi_divy_unit = {attribute,obj:'/nbi/divy', $
                     name:'units', $
                     data:'radians'}

    nbi_divz_desc = {attribute,obj:'/nbi/divz', $
                     name:'description', $
                     data:'Vertical divergences of the beam. One for each energy component'}
    nbi_divz_unit = {attribute,obj:'/nbi/divz', $
                     name:'units', $
                     data:'radians'}

    nbi_widy_desc = {attribute,obj:'/nbi/widy', $
                     name:'description', $
                     data:'Half width of the beam source grid'}
    nbi_widy_unit = {attribute,obj:'/nbi/widy', $
                     name:'units', $
                     data:'cm'}

    nbi_widz_desc = {attribute,obj:'/nbi/widz', $
                     name:'description', $
                     data:'Half height of the beam source grid'}
    nbi_widz_unit = {attribute,obj:'/nbi/widz', $
                     name:'units', $
                     data:'cm'}

    nbi_shape_desc = {attribute,obj:'/nbi/shape', $
                      name:'description', $
                      data:'Shape of the beam source grid: 1="rectangular", 2="circular"'}

    nbi_atts = [nbi_desc, nbi_cs, nbi_ds_desc,nbi_name_desc, $
                nbi_src_desc, nbi_src_unit, $
                nbi_axis_desc, nbi_axis_unit, $
                nbi_focy_desc, nbi_focy_unit, $
                nbi_focz_desc, nbi_focz_unit, $
                nbi_divy_desc, nbi_divy_unit, $
                nbi_divz_desc, nbi_divz_unit, $
                nbi_widy_desc, nbi_widy_unit, $
                nbi_widz_desc, nbi_widz_unit, $
                nbi_shape_desc]

    ;; Spectroscopic attributes
    spec_desc = {attribute,obj:'/spec', $
                 name:'description', $
                 data:'FIDA/BES Chord Geometry'}
    spec_cs = {attribute,obj:'/spec', $
               name:'coordinate_system', $
               data:'Right-handed cartesian'}

    spec_ds_desc = {attribute,obj:'/spec/data_source', $
                    name:'description', $
                    data:'Source of the chord geometry'}

    spec_nchan_desc = {attribute,obj:'/spec/nchan', $
                       name:'description', $
                       data:'Number of channels'}

    spec_system_desc = {attribute,obj:'/spec/system', $
                        name:'description', $
                        data:'Names of the different spectrocopic systems'}

    spec_lens_desc = {attribute,obj:'/spec/lens', $
                      name:'description', $
                      data:'Positions of the lenses'}
    spec_lens_unit = {attribute,obj:'/spec/lens', $
                      name:'units', $
                      data:'cm'}

    spec_axis_desc = {attribute,obj:'/spec/axis', $
                      name:'description', $
                      data:'Optical axis of the lines of sight: LOS(t) = lens + axis*t '}
    spec_axis_unit = {attribute,obj:'/spec/axis', $
                      name:'units', $
                      data:'cm'}

    spec_radius_desc = {attribute,obj:'/spec/radius', $
                        name:'description', $
                        data:'Line of sight radius at midplane or tangency point' }
    spec_radius_unit = {attribute,obj:'/spec/radius', $
                        name:'units', $
                        data:'cm'}

    spec_sigma_desc = {attribute,obj:'/spec/sigma_pi', $
                       name:'description', $
                       data:'Ratio of the intensities of the sigma and pi stark lines. Measured quantity'}

    spec_spot_desc = {attribute,obj:'/spec/spot_size', $
                      name:'description', $
                      data:'Radius of spot size'}
    spec_spot_unit = {attribute,obj:'/spec/spot_size', $
                      name:'units', $
                      data:'cm'}
 
    spec_atts = [spec_desc,spec_cs, spec_ds_desc, $
                 spec_nchan_desc, spec_system_desc, $
                 spec_lens_desc,spec_lens_unit, $
                 spec_axis_desc,spec_axis_unit, $
                 spec_radius_desc, spec_radius_unit, $
                 spec_sigma_desc, $
                 spec_spot_desc, spec_spot_unit]

    ;;NPA attributes
    npa_desc = {attribute,obj:'/npa', $
                name:'description', $
                data:'NPA Geometry'}

    npa_cs = {attribute,obj:'/npa', $
              name:'coordinate_system', $
              data:'Right-handed cartesian'}

    npa_ds_desc = {attribute,obj:'/npa/data_source', $
                   name:'description', $
                   data:'Source of the NPA geometry'}
    
    npa_nchan_desc = {attribute,obj:'/npa/nchan', $
                      name:'description', $
                      data:'Number of channels'}

    npa_system_desc = {attribute,obj:'/npa/system', $
                       name:'description', $
                       data:'Names of the different NPA systems'}


    npa_dshape_desc = {attribute,obj:'/npa/d_shape', $
                       name:'description', $
                       data:'Shape of the detector: 1="rectangular", 2="circular"'}

    npa_dcent_desc = {attribute,obj:'/npa/d_cent', $
                      name:'description', $
                      data:'Center of the detector'}
    npa_dcent_unit = {attribute,obj:'/npa/d_cent', $
                      name:'units', $
                      data:'cm'}

    npa_dtedge_desc = {attribute,obj:'/npa/d_tedge', $
                       name:'description', $
                       data:'Center of the detectors top edge'}
    npa_dtedge_unit = {attribute,obj:'/npa/d_tedge', $
                       name:'units', $
                       data:'cm'}

    npa_dredge_desc = {attribute,obj:'/npa/d_redge', $
                       name:'description', $
                       data:'Center of the detectors right edge'}
    npa_dredge_unit = {attribute,obj:'/npa/d_redge', $
                       name:'units', $
                       data:'cm'}

    npa_ashape_desc = {attribute,obj:'/npa/a_shape', $
                       name:'description', $
                       data:'Shape of the aperture: 1="rectangular", 2="circular"'}

    npa_acent_desc = {attribute,obj:'/npa/a_cent', $
                      name:'description', $
                      data:'Center of the aperture'}
    npa_acent_unit = {attribute,obj:'/npa/a_cent', $
                      name:'units', $
                      data:'cm'}

    npa_atedge_desc = {attribute,obj:'/npa/a_tedge', $
                       name:'description', $
                       data:'Center of the apertures top edge'}
    npa_atedge_unit = {attribute,obj:'/npa/a_tedge', $
                       name:'units', $
                       data:'cm'}

    npa_aredge_desc = {attribute,obj:'/npa/a_redge', $
                       name:'description', $
                       data:'Center of the apertures right edge'}
    npa_aredge_unit = {attribute,obj:'/npa/a_redge', $
                       name:'units', $
                       data:'cm'}

    npa_radius_desc = {attribute,obj:'/npa/radius', $
                       name:'description', $
                       data:'Line of sight radius at midplane or tangency point' }
    npa_radius_unit = {attribute,obj:'/npa/radius', $
                       name:'units', $
                       data:'cm'}

    npa_atts = [npa_desc, npa_cs, npa_ds_desc, $ 
                npa_nchan_desc, npa_system_desc, $
                npa_dshape_desc, npa_ashape_desc, $
                npa_dcent_desc, npa_dcent_unit, $
                npa_acent_desc, npa_acent_unit, $
                npa_dtedge_desc, npa_dtedge_unit, $
                npa_atedge_desc, npa_atedge_unit, $
                npa_dredge_desc, npa_dredge_unit, $
                npa_aredge_desc, npa_aredge_unit, $
                npa_radius_desc, npa_radius_unit ]

    atts = [root_atts, nbi_atts]
    geom = {nbi:nbi}

    if keyword_set(spec) then begin
        geom = create_struct(geom, "spec", spec)
        atts = [atts, spec_atts]
    endif

    if keyword_set(npa) then begin
        geom = create_struct(geom, "npa", npa)
        atts = [atts,npa_atts]
    endif

    write_hdf5, geom, filename=filename, atts=atts, /clobber

    if file_test(filename) then begin
        success, 'Geometry file created: '+filename
    endif else begin
        error, 'Geometry file creation failed.'
    endelse

END

PRO write_equilibrium, filename, plasma, fields
    ;+***
    ;+##`write_equilibrium, filename, plasma, fields`
    ;+Write MHD equilibrium values to a HDF5 file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the equilibrium file
    ;+
    ;+     **plasma**: Plasma structure
    ;+
    ;+     **fields**: Electromagnetic fields structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_equilibrium, filename, plasma, fields
    ;+```

    info, 'Writing equilibrium file...'

    root_atts = {attribute,obj:'/', $
                 name:'description', $
                 data:'Plasma Parameters and Electromagnetic Fields for FIDASIM'}

    ;; Plasma Attributes
    plasma_desc = {attribute,obj:'/plasma', $
                   name:'description', $
                   data:'Plasma Parameters'}
    
    plasma_cs = {attribute,obj:'/plasma', $
                 name:'coordinate_system', $
                 data:'Cylindrical'}

    plasma_ds_desc = {attribute, obj:'/plasma/data_source', $
                      name:'description', $
                      data:'Source of the plasma parameters'}

    plasma_time_desc = {attribute,obj:'/plasma/time', $
                        name:'description', $
                        data:'Time'}
    plasma_time_unit = {attribute,obj:'/plasma/time', $
                        name:'units', $
                        data:'s'}

    plasma_dene_desc = {attribute,obj:'/plasma/dene', $
                        name:'description', $
                        data:'Electron Number Density: Dene(r,z)'}
    plasma_dene_unit = {attribute,obj:'/plasma/dene', $
                        name:'units', $
                        data:'cm^-3'}

    plasma_te_desc = {attribute,obj:'/plasma/te', $
                      name:'description', $
                      data:'Electron Temperature: Te(r,z)'}
    plasma_te_unit = {attribute,obj:'/plasma/te', $
                      name:'units', $
                      data:'keV'}

    plasma_ti_desc = {attribute,obj:'/plasma/ti', $
                      name:'description', $
                      data:'Ion Temperature: Ti(r,z)'}
    plasma_ti_unit = {attribute,obj:'/plasma/ti', $
                      name:'units', $
                      data:'keV'}

    plasma_zeff_desc = {attribute,obj:'/plasma/zeff', $
                        name:'description', $
                        data:'Effective Nuclear Charge: Zeff(r,z)'}

    plasma_vr_desc = {attribute,obj:'/plasma/vr', $
                      name:'description', $
                      data:'Bulk plasma flow in the r-direction: Vr(r,z)'}
    plasma_vr_unit = {attribute,obj:'/plasma/vr', $
                      name:'units', $
                      data:'cm/s'}

    plasma_vt_desc = {attribute,obj:'/plasma/vt', $
                      name:'description', $
                      data:'Bulk plasma flow in the theta/torodial-direction: Vt(r,z)'}
    plasma_vt_unit = {attribute,obj:'/plasma/vt', $
                      name:'units', $
                      data:'cm/s'}

    plasma_vz_desc = {attribute,obj:'/plasma/vz', $
                      name:'description', $
                      data:'Bulk plasma flow in the z-direction: Vz(r,z)'}
    plasma_vz_unit = {attribute,obj:'/plasma/vz', $
                      name:'units', $
                      data:'cm/s'}

    plasma_nr_desc = {attribute,obj:'/plasma/nr', $
                      name:'description', $
                      data:'Number of R values'}
   
    plasma_nz_desc = {attribute,obj:'/plasma/nz', $
                      name:'description', $
                      data:'Number of Z values'}

    plasma_r_desc = {attribute,obj:'/plasma/r', $
                     name:'description', $
                     data:'Radius'}
    plasma_r_unit = {attribute,obj:'/plasma/r', $
                     name:'units', $
                     data:'cm'}

    plasma_z_desc = {attribute,obj:'/plasma/z', $
                     name:'description', $
                     data:'Z'}
    plasma_z_unit = {attribute,obj:'/plasma/z', $
                     name:'units', $
                     data:'cm'}

    plasma_r2d_desc = {attribute,obj:'/plasma/r2d', $
                       name:'description', $
                       data:'Radius grid: R(r,z)'}
    plasma_r2d_unit = {attribute,obj:'/plasma/r2d', $
                       name:'units', $
                       data:'cm'}

    plasma_z2d_desc = {attribute,obj:'/plasma/z2d', $
                       name:'description', $
                       data:'Z grid: Z(r,z)'}
    plasma_z2d_unit = {attribute,obj:'/plasma/z2d', $
                       name:'units', $
                       data:'cm'}

    plasma_mask_desc = {attribute,obj:'/plasma/mask',$
                         name:'description', $
                         data:'Boolean mask that indicates where' +$ 
                         ' the plasma parameters are well defined'}

    plasma_atts = [plasma_desc, plasma_cs, plasma_ds_desc, $
                   plasma_time_desc, plasma_time_unit, $
                   plasma_mask_desc, $
                   plasma_dene_desc, plasma_dene_unit, $
                   plasma_te_desc, plasma_te_unit, $
                   plasma_ti_desc, plasma_ti_unit, $
                   plasma_zeff_desc, $
                   plasma_vr_desc, plasma_vr_unit, $
                   plasma_vt_desc, plasma_vt_unit, $
                   plasma_vz_desc, plasma_vz_unit, $
                   plasma_nr_desc, plasma_nz_desc, $
                   plasma_r_desc, plasma_r_unit, $
                   plasma_z_desc, plasma_z_unit, $
                   plasma_r2d_desc, plasma_r2d_unit, $
                   plasma_z2d_desc, plasma_z2d_unit ]

    ;; Electromagnetic fields attributes
    fields_desc = {attribute,obj:'/fields', $
                   name:'description', $
                   data:'Electromagnetic Fields'}
    
    fields_cs = {attribute,obj:'/fields', $
                 name:'coordinate_system', $
                 data:'Cylindrical'}

    fields_ds_desc = {attribute,obj:'/fields/data_source', $
                      name:'description', $
                      data:'Source of the EM equilibrium'}

    fields_mask_desc = {attribute,obj:'/fields/mask',$
                        name:'description', $
                        data:'Boolean mask that indicates where' +$ 
                        ' the fields are well defined'}

    fields_time_desc = {attribute,obj:'/fields/time', $
                        name:'description', $
                        data:'Time'}
    fields_time_unit = {attribute,obj:'/fields/time', $
                        name:'units', $
                        data:'s'}

    fields_br_desc = {attribute,obj:'/fields/br', $
                      name:'description', $
                      data:'Magnetic field in the r-direction: Br(r,z)'}
    fields_br_unit = {attribute,obj:'/fields/br', $
                      name:'units', $
                      data:'T'}

    fields_bt_desc = {attribute,obj:'/fields/bt', $
                      name:'description', $
                      data:'Magnetic field in the theta/torodial-direction: Bt(r,z)'}
    fields_bt_unit = {attribute,obj:'/fields/bt', $
                      name:'units', $
                      data:'T'}

    fields_bz_desc = {attribute,obj:'/fields/bz', $
                      name:'description', $
                      data:'Magnetic field in the z-direction: Bz(r,z)'}
    fields_bz_unit = {attribute,obj:'/fields/bz', $
                      name:'units', $
                      data:'T'}

    fields_er_desc = {attribute,obj:'/fields/er', $
                      name:'description', $
                      data:'Electric field in the r-direction: Er(r,z)'}
    fields_er_unit = {attribute,obj:'/fields/er', $
                      name:'units', $
                      data:'V/m'}

    fields_et_desc = {attribute,obj:'/fields/et', $
                      name:'description', $
                      data:'Electric field in the theta/torodial-direction: Et(r,z)'}
    fields_et_unit = {attribute,obj:'/fields/et', $
                      name:'units', $
                      data:'V/m'}

    fields_ez_desc = {attribute,obj:'/fields/ez', $
                      name:'description', $
                      data:'Electric field in the z-direction: Ez(r,z)'}
    fields_ez_unit = {attribute,obj:'/fields/ez', $
                      name:'units', $
                      data:'V/m'}

    fields_nr_desc = {attribute,obj:'/fields/nr', $
                      name:'description', $
                      data:'Number of R values'}
   
    fields_nz_desc = {attribute,obj:'/fields/nz', $
                      name:'description', $
                      data:'Number of Z values'}

    fields_r_desc = {attribute,obj:'/fields/r', $
                     name:'description', $
                     data:'Radius'}
    fields_r_unit = {attribute,obj:'/fields/r', $
                     name:'units', $
                     data:'cm'}

    fields_z_desc = {attribute,obj:'/fields/z', $
                     name:'description', $
                     data:'Z'}
    fields_z_unit = {attribute,obj:'/fields/z', $
                     name:'units', $
                     data:'cm'}

    fields_r2d_desc = {attribute,obj:'/fields/r2d', $
                       name:'description', $
                       data:'Radius grid: R(r,z)'}
    fields_r2d_unit = {attribute,obj:'/fields/r2d', $
                       name:'units', $
                       data:'cm'}

    fields_z2d_desc = {attribute,obj:'/fields/z2d', $
                       name:'description', $
                       data:'Z grid: Z(r,z)'}
    fields_z2d_unit = {attribute,obj:'/fields/z2d', $
                       name:'units', $
                       data:'cm'}

    fields_atts = [fields_desc, fields_cs, fields_ds_desc, $
                   fields_mask_desc, $
                   fields_time_desc, fields_time_unit, $
                   fields_br_desc, fields_br_unit, $
                   fields_bt_desc, fields_bt_unit, $
                   fields_bz_desc, fields_bz_unit, $
                   fields_er_desc, fields_er_unit, $
                   fields_et_desc, fields_et_unit, $
                   fields_ez_desc, fields_ez_unit, $
                   fields_nr_desc, fields_nz_desc, $
                   fields_r_desc, fields_r_unit, $
                   fields_z_desc, fields_z_unit, $
                   fields_r2d_desc, fields_r2d_unit, $
                   fields_z2d_desc, fields_z2d_unit ]

    atts = [root_atts, plasma_atts, fields_atts]

    write_hdf5,["plasma","fields"],filename=filename, atts=atts, /clobber

    if file_test(filename) then begin
        success, 'Equilibrium file created: '+filename
    endif else begin
        error, 'Equilibrium file creation failed.'
    endelse
END

PRO write_distribution, filename, distri
    ;+***
    ;+##`write_distribution, filename, dist`
    ;+Write fast-ion distribution to a HDF5 file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the distribution file
    ;+
    ;+     **dist**: Fast-ion distribution structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_distribution, filename, distri
    ;+```

    info, 'Writing fast-ion distribution file...'

    root_atts = {attribute,obj:'/', $
                 name:'description', $
                 data:'Fast-ion distribution for FIDASIM'}

    cs_desc = {attribute,obj:'/', $
               name:'coordinate_system', $
               data:'Cylindrical'}

    ds_desc = {attribute,obj:'/data_source', $
               name:'description', $
               data:'Source of the fast-ion distribution'}    
    
    type_desc = {attribute,obj:'/type', $
                 name:'description', $
                 data:'Distribution type: '+ $
                      '1="Guiding Center Density Function", '+ $
                      '2="Guiding Center Monte Carlo", '+ $
                      '3="Full Orbit Monte Carlo"'}

    time_desc = {attribute,obj:'/time', $
                 name:'description', $
                 data:'Distribution time'}
    time_unit = {attribute,obj:'/time', $
                 name:'units', $
                 data:'s'}

    type = distri.type

    if type eq 1 then begin
        nen_desc = {attribute,obj:'/nenergy',$
                    name:'description', $
                    data:'Number of energy values'}
        np_desc = {attribute,obj:'/npitch', $
                   name:'description', $
                   data:'Number of pitch values'}

        energy_desc = {attribute,obj:'/energy', $
                       name:'description',$
                       data:'Energy'}
        energy_unit = {attribute,obj:'/energy', $
                       name:'units', $
                       data:'keV'}
 
        pitch_desc = {attribute,obj:'/pitch', $
                       name:'description',$
                       data:'Pitch: p = v_parallel/v  w.r.t. the magnetic field'}

        f_desc = {attribute,obj:'/f', $
                  name:'description', $
                  data:'Fast-ion density function: F(E,p,R,Z)'}
        f_unit = {attribute,obj:'/f', $
                  name:'units', $
                  data:'fast-ions/(dE*dP*cm^3)'} 

        denf_desc = {attribute,obj:'/denf', $
                     name:'description', $
                     data:'Fast-ion density: Denf(r,z)'}
        denf_unit = {attribute,obj:'/denf', $
                     name:'units', $
                     data:'cm^-3'} 

        nr_desc = {attribute,obj:'/nr', $
                   name:'description', $
                   data:'Number of R values'}
   
        nz_desc = {attribute,obj:'/nz', $
                   name:'description', $
                   data:'Number of Z values'}

        r_desc = {attribute,obj:'/r', $
                  name:'description', $
                  data:'Radius'}
        r_unit = {attribute,obj:'/r', $
                  name:'units', $
                  data:'cm'}

        z_desc = {attribute,obj:'/z', $
                  name:'description', $
                  data:'Z'}
        z_unit = {attribute,obj:'/z', $
                  name:'units', $
                  data:'cm'}

        r2d_desc = {attribute,obj:'/r2d', $
                    name:'description', $
                    data:'Radius grid: R(r,z)'}
        r2d_unit = {attribute,obj:'/r2d', $
                    name:'units', $
                    data:'cm'}

        z2d_desc = {attribute,obj:'/z2d', $
                    name:'description', $
                    data:'Z grid: Z(r,z)'}
        z2d_unit = {attribute,obj:'/z2d', $
                    name:'units', $
                    data:'cm'}

        atts = [root_atts, cs_desc, ds_desc, $
                type_desc, time_desc,time_unit, $
                nen_desc, np_desc, $
                energy_desc, energy_unit, $
                pitch_desc, $
                f_desc, f_unit, $
                denf_desc, denf_unit, $
                nr_desc, nz_desc, $
                r_desc, r_unit, $
                z_desc, z_unit, $
                r2d_desc, r2d_unit, $
                z2d_desc, z2d_unit]
        
    endif else begin

        np_desc = {attribute,obj:'/nparticle', $
                   name:'description', $
                   data:'Number of MC particles'}

        nc_desc = {attribute,obj:'/nclass', $
                   name:'description', $
                   data:'Number of orbit classes'}

        r_desc = {attribute,obj:'/r', $
                  name:'description', $
                  data:'R position of a MC particle'}
        r_unit = {attribute,obj:'/r', $
                  name:'units', $
                  data:'cm'}

        z_desc = {attribute,obj:'/z', $
                  name:'description', $
                  data:'Z position of a MC particle'}
        z_unit = {attribute,obj:'/z', $
                  name:'units', $
                  data:'cm'}

        w_desc = {attribute,obj:'/weight', $
                  name:'description', $
                  data:'Weight of a MC particle: sum(weight) = # of fast-ions '}
        w_unit = {attribute,obj:'/weight', $
                  name:'units', $
                  data:'fast-ions/particle'}

        c_desc = {attribute,obj:'/class', $
                  name:'description', $
                  data:'Orbit class of a MC particle: class in Set(1:nclass)'}

        if type eq 2 then begin
            energy_desc = {attribute,obj:'/energy', $
                           name:'description', $
                           data:'Energy of a MC particle'}
            energy_unit = {attribute,obj:'/energy', $
                           name:'units', $
                           data:'keV'}

            pitch_desc = {attribute,obj:'/pitch', $
                          name:'description', $
                          data:'Pitch of a MC particle: p = v_parallel/v  w.r.t. the magnetic field'}
            type_atts = [energy_desc,energy_unit,pitch_desc]
        endif else begin
            vr_desc = {attribute,obj:'/vr', $
                       name:'description', $
                       data:'Radial velocity of a MC particle'}
            vr_unit = {attribute,obj:'/vr', $
                       name:'units', $
                       data:'cm/s'}
            vt_desc = {attribute,obj:'/vt', $
                       name:'description', $
                       data:'Torodial velocity of a MC particle'}
            vt_unit = {attribute,obj:'/vt', $
                       name:'units', $
                       data:'cm/s'}
            vz_desc = {attribute,obj:'/vz', $
                       name:'description', $
                       data:'Z velocity of a MC particle'}
            vz_unit = {attribute,obj:'/vz', $
                       name:'units', $
                       data:'cm/s'}
            type_atts = [vr_desc,vr_unit,vt_desc,vt_unit,vz_desc,vz_unit]
        endelse

        atts = [root_atts, cs_desc, ds_desc, $
                type_desc, time_desc,time_unit, $
                np_desc, nc_desc, $
                r_desc, r_unit, $
                z_desc, z_unit, $
                w_desc, w_unit, $
                c_desc, type_atts] 
    endelse

    write_hdf5, distri, filename=filename,atts=atts, /clobber
   
    if file_test(filename) then begin
        success, 'Distribution file created: '+filename
    endif else begin
        error, 'Distribution file creation failed.'
    endelse
END

PRO prefida,inputs,grid,nbi,plasma,fields,fbm,spec=spec,npa=npa
    ;+***
    ;+##`prefida, inputs, grid, nbi, plasma, fields, dist, spec=spec, npa=npa`
    ;+Checks FIDASIM inputs and writes FIDASIM input HDF5 files
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
    check_dist, inputs, grid, fbm, err_status
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
