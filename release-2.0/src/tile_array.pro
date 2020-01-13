FUNCTION tile_array, arr, ncol, nrow
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
