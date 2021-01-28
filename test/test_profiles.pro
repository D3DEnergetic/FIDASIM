FUNCTION test_profiles,filename,grid,rhogrid

    prof = read_ncdf(filename)

    impurity_charge=6
    nthermal = 1
    species_mass = [2.01410178d0]

    rho = double(prof.rho)
    dene = interpol(prof.dene*1.0d-6,rho,rhogrid) ;;cm^-3
    ti = interpol(prof.ti*1.0d-3,rho,rhogrid) ;;keV
    te = interpol(prof.te*1.0d-3,rho,rhogrid) ;;keV
    zeff = dene*0.d0 + 1.5
    denimp = dene*(zeff - 1)/(impurity_charge*(impurity_charge-1))
    deni = dene - impurity_charge*denimp
    vt = grid.r2d*interpol(prof.omega*1.0d0,rho,rhogrid) ;;cm/s
    vr = 0.d0*vt ;;cm/s
    vz = 0.d0*vt ;;cm/s
    denn = zeff*0 + 1.0d8
    max_rho = max(abs(rho))

    s = size(rhogrid,/dim)
    mask = intarr(s[0],s[1])
    w=where(rhogrid le max_rho) ;where we have profiles
    mask[w] = 1

    nr = n_elements(dene[*,0])
    nz = n_elements(dene[0,*])
    dene=rebin(dene,nr,nz,grid.nphi)
    denimp=rebin(denimp,nr,nz,grid.nphi)
    deni=reform(deni,nthermal,nr,nz,grid.nphi)
    denn=rebin(denn,nr,nz,grid.nphi)
    mask=rebin(mask,nr,nz,grid.nphi)
    te=rebin(te,nr,nz,grid.nphi)
    ti=rebin(ti,nr,nz,grid.nphi)
    vt=rebin(vt,nr,nz,grid.nphi)
    vr=rebin(vr,nr,nz,grid.nphi)
    vz=rebin(vz,nr,nz,grid.nphi)
    zeff=rebin(zeff,nr,nz,grid.nphi)

    ;;SAVE IN PROFILE STRUCTURE
    profiles={time:1.d0,data_source:filename,mask:mask, nthermal:nthermal,$
              denimp:denimp, deni:deni,species_mass:species_mass,impurity_charge:impurity_charge, $
              te:te,ti:ti,vr:vr,vt:vt,vz:vz,dene:dene,zeff:zeff,denn:denn}

    return,profiles
END
