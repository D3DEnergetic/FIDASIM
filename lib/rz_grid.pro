FUNCTION rz_grid,rmin,rmax,nr,zmin,zmax,nz
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
    ;+IDL> grid = rz_grid(0,200.0,200,-100,100,200)
    ;+```

    dr = (rmax-rmin)/(nr-1)
    dz = (zmax-zmin)/(nz-1)
    r = rmin + dr*dindgen(nr)
    z = zmin + dz*dindgen(nz)

    r2d = r # replicate(1,nz)
    z2d = replicate(1,nr) # z


    grid = {r2d:r2d,z2d:z2d,r:r,z:z,nr:nr,nz:nz}
    
    return, grid
END
