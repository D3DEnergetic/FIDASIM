FUNCTION test_profiles,filename,grid,flux

    ;; profiles structure
    ;;** Structure <83e8518>, 7 tags, length=5768, data length=5764, refs=1:
    ;;   DATA_SOURCE     STRING
    ;;   TIME            DOUBLE                 [s]
    ;;   MASK            INT
    ;;   VR              DOUBLE    Array[nr,nz] [cm/s]
    ;;   VT              DOUBLE    Array[nr,nz] [cm/s]
    ;;   VZ              DOUBLE    Array[nr,nz] [cm/s]
    ;;   TI              DOUBLE    Array[nr,nz] [keV]
    ;;   TE              DOUBLE    Array[nr,nz] [keV]
    ;;   DENE            DOUBLE    Array[nr,nz] [cm^-3]
    ;;   ZEFF            DOUBLE    Array[nr,nz]
    ;;   DENN            DOUBLE    Array[nr,nz] [cm^-3]


    prof = read_ncdf(filename)

    rho = double(prof.rho)
    dene = interpol(prof.dene*1.0d-6,rho,flux) ;;cm^-3
    ti = interpol(prof.ti*1.0d-3,rho,flux) ;;keV
    te = interpol(prof.te*1.0d-3,rho,flux) ;;keV
    zeff = dene*0.d0 + 1.5
    vt = grid.r2d*interpol(prof.omega*1.0d0,rho,flux) ;;cm/s
    vr = 0.d0*vt ;;cm/s
    vz = 0.d0*vt ;;cm/s
    denn = zeff*0 + 1.0d8
    max_flux = max(abs(rho))

    s = size(flux,/dim)
    mask = intarr(s[0],s[1])
    w=where(flux le max_flux) ;where we have profiles
    mask[w] = 1

    nr = n_elements(dene[*,0])
    nz = n_elements(dene[0,*])
    dene=rebin(dene,nr,nz,grid.nphi)
    denn=rebin(denn,nr,nz,grid.nphi)
    mask=rebin(mask,nr,nz,grid.nphi)
    te=rebin(te,nr,nz,grid.nphi)
    ti=rebin(ti,nr,nz,grid.nphi)
    vt=rebin(vt,nr,nz,grid.nphi)
    vr=rebin(vr,nr,nz,grid.nphi)
    vz=rebin(vz,nr,nz,grid.nphi)
    zeff=rebin(zeff,nr,nz,grid.nphi)

    ;;SAVE IN PROFILE STRUCTURE
	profiles={time:1.d0,data_source:filename,mask:mask, $
              te:te,ti:ti,vr:vr,vt:vt,vz:vz,dene:dene,zeff:zeff,denn:denn}

	return,profiles
END
