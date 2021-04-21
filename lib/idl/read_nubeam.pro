FUNCTION grid_fbm,r2d,z2d,fbm,fdens,rout,zout
    compile_opt idl2, logical_predicate

    s = size(fbm,/dim)
    nenergy = s[0]
    npitch = s[1]
    nr = n_elements(rout)
    nz = n_elements(zout)
    nr2d = n_elements(r2d)
    dim = [nr,nz]
    delta = abs([rout[1]-rout[0],zout[1]-zout[0]])
    start = [min(rout),min(zout)]

    triangulate, r2d, z2d, tr

    r2dt = r2d[tr[*]]
    z2dt = z2d[tr[*]]
    linTr = lindgen(size(tr,/dim))
    index = lindgen(n_elements(tr))/3*3
    tr_num = round(griddata(r2dt,z2dt,float(index),triangles=linTr,/linear,$
             start=start,delta=delta,dimension=dim))

    wts = ptrarr(3)
    for i=0,2 do begin
        w = griddata(r2dt,z2dt,(lindgen(n_elements(r2dt)) mod 3) eq i,triangles=linTr,/linear,$
            start=start,delta=delta,dimension=dim)
        wts[i] = ptr_new(w,/no_copy)
    endfor

    denf = dblarr(nr,nz)
    for i=0,2 do begin
        denf = denf + fdens[tr[tr_num+i]]*(*wts[i])
    endfor
    denf = denf > 0

    fbm_grid = dblarr(nenergy,npitch,nr,nz)
    for i=0,nenergy-1 do begin
        for j=0,npitch-1 do begin
            for k=0,2 do begin
                fbm_grid[i,j,*,*] = fbm_grid[i,j,*,*] + fbm[i,j,tr[tr_num+k]]*(*wts[k])
            endfor
        endfor
    endfor
    fbm_grid = fbm_grid > 0

    ;; Catch points outside of triangulation
    tr_ind = fix(griddata(r2d,z2d,indgen(nr2d),triangles=tr,/nearest_neighbor, $
                      start=start,delta=delta,dimension=dim))
    denf_nn = fdens[tr_ind]
    fbm_grid_nn = fbm[*,*,tr_ind]

    w = where(denf le 0.0,nw)
    if nw ne 0 then begin
        denf[w] = denf_nn[w]
        inds=array_indices(denf,w)
        for i=0,nw-1 do begin
            fbm_grid[*,*,inds[0,i],inds[1,i]]=fbm_grid_nn[*,*,w[i]]
        endfor
    end

    return, {denf:denf, fbm:fbm_grid}
END

FUNCTION read_nubeam,filename,grid,btipsign=btipsign,e_range=e_range,p_range=p_range,species=species
    ;+#read_nubeam
    ;+Reads NUBEAM fast-ion distribution function
    ;+***
    ;+##Arguments
    ;+    **filename**: NUBEAM guiding center fast-ion distribution function file e.g. 159245H01_fi_1.cdf
    ;+
    ;+    **grid**: Interpolation grid
    ;+
    ;+##Keyword Arguments
    ;+    **btipsign**: Sign of the dot product of the magnetic field and plasma current
    ;+
    ;+    **e_range**: Energy range to consider
    ;+
    ;+    **p_range**: Pitch range to consider
    ;+
    ;+    **species**: Fast-ion species number. Defaults to 1
    ;+
    ;+##Return Value
    ;+Distribution structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> dist = read_nubeam("./159245H02_fi_1.cdf",grid,btipsign=-1)
    ;+```

    if not keyword_set(btipsign) then btipsign = -1

    cdftest=findfile(filename)

    if cdftest[0] eq '' then begin
        err=1
        goto,GET_OUT
    endif

    if not keyword_set(species) then species = 1

    sstr = string((read_ncdf(filename,vars=["SPECIES_"+strcompress(species,/r)])).(1))
    print, 'Species: ' + sstr
    vars = read_ncdf(filename,vars=["TIME","R2D","Z2D","E_"+sstr,"A_"+sstr, $
                                    "F_"+sstr,"RSURF","ZSURF","BMVOL"])
    ngrid=n_elements(vars.r2d)

    ;;-------------Convert eV-> keV
    time = vars.time
    r2d = vars.r2d
    z2d = vars.z2d
    rsurf = vars.rsurf
    zsurf = vars.zsurf
    sstr = strupcase(sstr)
    index = where(tag_names(vars) EQ "A_"+sstr)
    pitch = vars.(index[0])
    index = where(tag_names(vars) EQ "E_"+sstr)
    energy=vars.(index[0])*1.0d-3          ;; fidasim needs energy in kev
    index = where(tag_names(vars) EQ "F_"+sstr)
    fbm=vars.(index[0])*1.0d3              ;; now, this needs to be corrected
    ;; as we now calculate with fast-ions/omega/keV/cm^3
    ;;------------Convert d_omega --> pitch
    ;; Fast-ion distribution is given as a function of cm^3, energy
    ;; and d_omega/4Pi. omega is the solild angle in 3D velocity space. In
    ;; order to transform this to a function depending on pitch instead
    ;; of d_omega/4PI, one has to multiply by 0.5!
    fbm*=0.5
    ;; make sure that fbm is >=0:
    fbm>=0.
    ;;loading finished

    ;; TRANSP defines the pitch along the current direction. In
    ;; contrast, FIDASIM uses pitch along the B-field! Therefore,
    ;; reverse the pitch coordinate in fbm if B and J are anti-parallel!
    if btipsign lt 0 then begin
        npitch=n_elements(pitch)
        index=npitch-(indgen(npitch)+1)
        fbm[*,*,*]=fbm[*,index,*]
    endif

    if not keyword_set(e_range) then begin
        e_range = [min(energy),max(energy)]
    endif

    if not keyword_set(p_range) then begin
        p_range = [min(pitch),max(pitch)]
    endif

    ;;----------- select energy range -------
    index=where(energy ge e_range[0] and energy le e_range[1],nenergy)
    energy=energy[index]
    fbm=fbm[index,*,*]
    dE      = energy[1] - energy[0]
    emin=(float(energy[0])         - float(0.5*dE))>0.
    emax=float(energy[nenergy-1]) + float(0.5*dE)
    print, 'Energy min/max:', emin,emax

    ;; --------- select Pitch range --------
    index=where(pitch ge p_range[0] and pitch le p_range[1],npitch)
    pitch=pitch[index]
    fbm=fbm[*,index,*]
    dP  = abs(pitch[1]  - pitch[0])
    pmin=(float(pitch[0])       - float(0.5*dP))>(-1)
    pmax=(float(pitch[npitch-1])+ float(0.5*dP))<1
    print, 'Pitch  min/max:', pmin,pmax

    ;; ------map fdens on FIDASIM grid and sort out
    ;; ------points outside the separatrix
    nr=grid.nr
    nz=grid.nz
    nphi=grid.nphi
    rgrid=grid.r
    zgrid=grid.z
    dr = abs(rgrid[1]-rgrid[0])
    dz = abs(zgrid[1]-zgrid[0])

    ;; FBM & DENF
    fdens=total(reform(total(fbm,1)),1)*dE*dP
    ntot = total(fdens*vars.bmvol)
    print, 'Ntotal in phase space: ',ntot
    fstr = grid_fbm(r2d,z2d,fbm,fdens,rgrid,zgrid)
    denf = fstr.denf
    fbm_grid=fstr.fbm

    ;; sort out positions more than 2 cm outside the separatrix
    rmaxis=mean(rsurf[*,0])
    zmaxis=mean(zsurf[*,0])
    rsep=rsurf[*,(size(rsurf))[2]-1]
    zsep=zsurf[*,(size(rsurf))[2]-1]
    x_bdry = rsep - rmaxis
    y_bdry = zsep - zmaxis
    r_bdry = sqrt(x_bdry^2 + y_bdry^2)
    theta  = atan(y_bdry,x_bdry)
    ;; -- sort and remove identical values --
    index = uniq(theta,sort(theta))
    theta           = theta[index]
    r               = r_bdry[index]
    ;; --- make theta periodic
    n               = n_elements(r)
    r_bdry          = fltarr(n+2)
    theta_bdry      = fltarr(n+2)
    theta_bdry[1:n] = theta
    r_bdry[1:n]     = r
    r_bdry[0]       = r[n-1]
    theta_bdry[0]   = theta[n-1] - 2.*!pi
    r_bdry[n+1]     = r[0]
    theta_bdry[n+1] = theta[0] + 2.*!pi
    ;; -- express (r_pts,z_pts) in (r,theta) coordinates --
    x_pts = grid.r2d - rmaxis
    y_pts = grid.z2d - zmaxis
    r_pts = sqrt(x_pts^2 + y_pts^2)
    theta_pts=atan(y_pts,x_pts)
    ;; -- interpolate to get the radial position of the boundary
    ;;    evaluated at theta = theta_pts --
    index=sort(theta_pts)
    mapped = interpol(r_bdry,theta_bdry,theta_pts[index])
    r_boundary = theta_pts*0.d
    r_boundary[index]=mapped
    index = where(r_pts gt r_boundary+2., nind)
    if nind gt 0 then begin
        indices=array_indices(r_pts,index)
        for i=0,(size(indices))[2]-1 do begin
            fbm_grid[*,*,indices[0,i],indices[1,i]]=0.
            denf[indices[0,i],indices[1,i]]=0.
        endfor
    endif

    ;; enforce correct normalization
    ntot_denf = (2*!dpi*dr*dz)*total(rgrid*total(denf,2))
    denf = denf*(ntot/ntot_denf)
    ntot_fbm = (2*!dpi*dr*dz*dE*dP)*total(rgrid*total(total(total(fbm_grid,1),1),2))
    fbm_grid = fbm_grid*(ntot/ntot_fbm)

    denf=rebin(denf,nr,nz,nphi)
    fbm_grid=rebin(fbm_grid,nenergy,npitch,nr,nz,nphi)

    fbm_struct={type:1,time:time,nenergy:fix(nenergy),energy:energy,npitch:fix(npitch),$
                pitch:pitch,f:fbm_grid,denf:denf,data_source:file_expand_path(filename)}

    return, fbm_struct
    GET_OUT:
END
