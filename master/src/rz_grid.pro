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
    ;+IDL> grid = rz_grid(100.d0,240.d0, 70, -100.d0,100.d0, 100, phimin=4*!dpi/3, phimax=5*!dpi/3, nphi=5)
    ;+```

    if not keyword_set(phimin) then phimin = 0.0 ;rad
    if not keyword_set(phimax) then phimax = 0.0 ;rad
    if not keyword_set(nphi) then nphi = 1

    dr = double(rmax-rmin)/nr
    dz = double(zmax-zmin)/nz
    dphi = double(phimax-phimin)/(nphi)
    r = rmin + dr*dindgen(nr)
    z = zmin + dz*dindgen(nz)
    phi = phimin + dphi*dindgen(nphi)

    r2d = r # replicate(1,nz)
    z2d = replicate(1,nr) # z


    grid = {r2d:r2d,z2d:z2d,r:r,z:z,phi:phi,nr:nr,nz:nz,nphi:nphi}
    
    return, grid
END
