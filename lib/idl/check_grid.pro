PRO check_grid, grid
    ;+#check_grid
    ;+Checks if interpolation grid structure is valid
    ;+***
    ;+##Input Arguments
    ;+     **grid**: Interpolation grid structure
    ;+ 
    ;+##Example Usage
    ;+```idl
    ;+IDL> check_grid, grid
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

    w = where("nphi" eq strlowcase(TAG_NAMES(grid)),nw)
    if nw eq 0 then begin
        info,'"nphi" is missing from the interpolation grid, assuming axisymmetry'
        nphi = 1
    endif else begin
        nphi = grid.nphi
    endelse

    nr = grid.nr
    nz = grid.nz
    zero_int = {dims:0,type:'INT'}
    schema = {nr:zero_int, nz:zero_int, nphi:zero_int, $
              r2d:{dims:[nr,nz], type:'DOUBLE'}, $
              z2d:{dims:[nr,nz], type:'DOUBLE'}, $
              r:{dims:[nr], type:'DOUBLE'}, $
              z:{dims:[nz], type:'DOUBLE'}, $
              phi:{dims:[nphi], type:'DOUBLE'} }

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

    w = where((indgen(nphi) eq sort(phi)) ne 1, nw)
    if nw ne 0 then begin
        error,'phi is not in ascending order'
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
        error,'Invalid interpolation grid. Exiting...',/halt
    endif else begin
        success,'Interpolation grid is valid'
    endelse

END
