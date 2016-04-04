FUNCTION rz_grid,rmin,rmax,nr,zmin,zmax,nz

    dr = (rmax-rmin)/(nr-1)
    dz = (zmax-zmin)/(nz-1)
    r = rmin + dr*dindgen(nr)
    z = zmin + dz*dindgen(nz)

    r2d = r # replicate(1,nz)
    z2d = replicate(1,nr) # z


    grid = {r2d:r2d,z2d:z2d,r:r,z:z,nr:nr,nz:nz}
    
    return, grid
END
