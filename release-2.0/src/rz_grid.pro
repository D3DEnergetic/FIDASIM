FUNCTION rz_grid,rmin,rmax,nr,zmin,zmax,nz,phimin=phimin,phimax=phimax,nphi=nphi
    ;+#rz_grid
    ;+Creates interpolation grid
    ;+***
    ;+##Arguments
    ;+    **rmin**: Minimum radius [cm]
    ;+
    ;+    **rmax**: Maximum radius [cm]
    ;+
    ;+    **nr**: Number of radii
    ;+
    ;+    **zmin**: Minimum Z value [cm]
    ;+
    ;+    **zmax**: Maximum Z value [cm]
    ;+
    ;+    **nz**: Number of z values
    ;+
    ;+##Return Value
    ;+Interpolation grid structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> grid = rz_grid(100.d0,240.d0,70,-100.d0,100.d0,100,phimin=-!dpi/6,phimax=!dpi/6,nphi=5)
    ;+```

    if not keyword_set(phimin) then phimin = 0.0 ;rad
    if not keyword_set(phimax) then phimax = 0.0 ;rad
    if not keyword_set(nphi) then nphi = 1

    dr = double(abs(rmax-rmin))/(nr-1)
    dz = double(abs(zmax-zmin))/(nz-1)
    if nphi eq 1 then begin
        dphi = 0.0
    endif else begin
        dphi = double(abs(phimax-phimin))/(nphi-1)
    endelse

    r = rmin + dr*dindgen(nr)
    z = zmin + dz*dindgen(nz)
    phi = phimin + dphi*dindgen(nphi)

    r2d = r # replicate(1,nz)
    z2d = replicate(1,nr) # z

    grid = {r2d:r2d,z2d:z2d,r:r,z:z,phi:phi,nr:nr,nz:nz,nphi:nphi}

    return, grid
END
