FUNCTION nstx_equil,inputs,grid,det

    equil={err:1}
    ;; Get eqdsk
    time_str='00000'+strtrim(string(long(inputs.time*1000)),1)
    time_str=strmid(time_str,4,/reverse_offset)
    shot_str=strtrim(string(inputs.shot),1)
    profile_str='g'+shot_str+'.'+time_str
    slash=strmid(inputs.profile_dir,0,1,/reverse_offset)
    if slash ne '/' then begin
       gfile=inputs.profile_dir+'/'+profile_str
    endif else begin
       gfile=inputs.profile_dir+profile_str
    endelse    
    gfiletest=findfile(gfile)
    if gfiletest ne '' then begin
        print,'Restoring equilibrium from gfile'
        restore,gfile
    endif else begin
        print,'gfile does not exist'
        goto, GET_OUT
    endelse

    rhogrid=rho_rz(g,grid.r_grid/100.,grid.w_grid/100.,/do_linear)

    calculate_bfield,bp,br,bphi,bz1,g

    ;; Get radial electric field on efit's grid from potential
    ;; epoten is on a grid of equally spaced points in psi from g.ssimag to g.ssibry
    dpsi=(g.ssibry-g.ssimag)/(n_elements(g.epoten)-1)
    psi=g.ssimag + dpsi*findgen(n_elements(g.epoten))
    npot=n_elements(g.epoten)
    epot=replicate(g.epoten[npot-1],n_elements(g.r),n_elements(g.z))
    for i=0l,n_elements(g.r)-1 do begin 
        for j=0l,n_elements(g.z)-1 do begin
            psi1=g.psirz[i,j]
            dum=min(abs(psi1-psi),kpsi)
            if kpsi ne npot-1 then epot[i,j]=spline(psi,g.epoten,[psi1])
        endfor
    endfor
    ; E = - grad(Phi)    EFIT units should be V/m
    er=-(shift(epot,-1,0) - shift(epot,1,0))/(g.r[2]-g.r[0])
    ez1=-(shift(epot,0,-1) - shift(epot,0,1))/(g.z[2]-g.z[0])
    
    ;; Interpolate cylindrical fields onto (x,y,z) mesh then rotate vectors
    bx=dblarr(grid.nx,grid.ny,grid.nz) & by=bx & bz=bx
    ex=dblarr(grid.nx,grid.ny,grid.nz) & ey=ex & ez=ex  
    
    for i=0L,grid.nx-1 do for j=0L,grid.ny-1 do for k=0L,grid.nz-1 do begin
        rgrid=(.01*grid.r_grid[i,j,k] - g.r[0])/(g.r[1]-g.r[0]) ; in grid units
        zgrid=(.01*grid.w_grid[i,j,k] - g.z[0])/(g.z[1]-g.z[0])    ; WWH 3/31/07
        bcylr=interpolate(br,[rgrid],[zgrid])
        ecylr=interpolate(er,[rgrid],[zgrid])
        bcylphi=interpolate(bphi,[rgrid],[zgrid])
        ez[i,j,k]=interpolate(ez1,[rgrid],[zgrid])
        bz[i,j,k]=interpolate(bz1,[rgrid],[zgrid])
        cph=cos(grid.phi_grid[i,j,k]) & sph=sin(grid.phi_grid[i,j,k])
        bx[i,j,k]=(cph*bcylr - sph*bcylphi)
        by[i,j,k]=(sph*bcylr + cph*bcylphi)
        ex[i,j,k]=cph*ecylr
        ey[i,j,k]=sph*ecylr
    endfor
    
    if inputs.calc_brems eq 0 then begin
        ;;GET RHO VALUES ALONG LINE OF SIGHT
        ds=.3    ; step size (cm)
        ns=4000  ; maximum number of steps

        ;; Sightlines
        nchan=det.nchan
        x0=det.xlens   ; cm
        y0=det.ylens
        z0=det.zlens
        v=fltarr(3,nchan)
        for i=0,nchan-1 do begin
            v[0,i]=det.xlos[i]-x0[i]
            v[1,i]=det.ylos[i]-y0[i]
            v[2,i]=det.zlos[i]-z0[i]
            v[*,i]=v[*,i]/sqrt(v[0,i]^2+v[1,i]^2+v[2,i]^2)
        endfor
        rhospath=dblarr(ns,nchan)
        for i=0,nchan-1 do begin
            ;equations for lines in x,y,z are
            x=x0[i] + v[0,i]*ds*findgen(ns)  ; cm
            y=y0[i] + v[1,i]*ds*findgen(ns)
            z=z0[i] + v[2,i]*ds*findgen(ns)
    
            rmajvals=(x^2.+y^2.)^.5  ; major radii of all points along ray
            rhospath[*,i]=rho_rz(g,rmajvals/100.,z/100.,/do_linear)
        endfor  ; channel loop

        rho_chords={rhos:rhospath,ds:ds}
    endif else rho_chords={rhos:0,ds:0}

    equil={g:g,rho_grid:rhogrid,rho_chords:rho_chords,bx:bx,by:by,bz:bz,ex:ex,ey:ey,ez:ez,err:0}
    GET_OUT:
    return,equil
END
