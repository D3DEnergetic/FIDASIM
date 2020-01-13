FUNCTION is_left, p0, p1, p2
    return, ((p1[0] - p0[0])*(p2[1] - p0[1]) - (p2[0] - p0[0])*(p1[1] - p0[1]))
END

FUNCTION in_vessel, rwall, zwall, r, z
    ;+#in_vessel
    ;+Calculates whether a given r,z is within the vessel
    ;+***
    ;+##Arguments
    ;+    **rwall**: r values of the wall
    ;+
    ;+    **zwall**: z values of the wall
    ;+
    ;+    **r**: r values
    ;+
    ;+    **z**: z values
    ;+
    ;+##Return Value
    ;+    0 for false, 1 for true
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> inv = in_vessel(rwall,zwall,r,z)
    ;+```

    if n_elements(rwall) ne n_elements(zwall) then begin
        error, "rwall and zwall must have the same number of elements"
        stop
    endif
    nv = n_elements(rwall)
    v = [[reform(rwall,nv)],[reform(zwall,nv)]]

    if n_elements(r) ne n_elements(z) then begin
        error, "r and z must have the same number of elements"
        stop
    endif
    np = n_elements(r)
    dims = size(r,/d)
    rr = reform(r,np)
    zz = reform(z,np)

    inv = intarr(np)
    for ip=0, np-1 do begin
        p = [rr[ip], zz[ip]]

        wn = 0
        for i = 0, nv-1 do begin
            i1 = (i + 1) mod nv
            if v[i,1] le p[1] then begin
                if v[i1,1] gt p[1] then begin
                    if is_left(v[i,*],v[i1,*],p) gt 0.0 then begin
                        wn = wn + 1
                    endif
                endif
            endif else begin
                if v[i1,1] le p[1] then begin
                    if is_left(v[i,*],v[i1,*],p) lt 0.0 then begin
                        wn = wn - 1
                    endif
                endif
            endelse
        endfor

        inv[ip] = wn ne 0
    endfor

    return, reform(inv,dims)
END
