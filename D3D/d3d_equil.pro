FUNCTION sinterpol,v,x,u,sortt=sortt,_extra=_extra
    if n_elements(sortt) lt 1 then sortt=0

    if sortt then begin
        ind=sort(X)
    endif else begin
        ind=lindgen(n_elements(x))
    endelse


    return,interpol(v[ind],x[ind],u,_extra=_extra)
END

FUNCTION interpolatexy,x,y,z,xu,yu,_extra=_extra

	;takes in vectors of x and y coordinates for z(x,y) data set
	;interpolates z onto (xu,yu) coordinates has same keywords as 
	;interpolate


	indx=sinterpol(indgen(n_elements(x)),x,xu,sortt=1)
	indy=sinterpol(indgen(n_elements(y)),y,yu,sortt=1)

	int=interpolate(z,indx,indy,_extra=_extra)

	return,int

END

FUNCTION d3d_equil,inputs,grid,det

	equil={err:1}
	;; Get eqdsk
	if inputs.gfile ne '' then begin
		gfiletest=findfile(inputs.gfile)
		if gfiletest[0] eq '' then begin 
			print,'FATAL ERROR in D3D_EQUIL: gfile ', inputs.gfile, ' not found'
			goto,GET_OUT
		endif
		g=readg(inputs.gfile) 
	endif else begin
		g=readg(inputs.shot,inputs.time*1000,RUNID=inputs.equil,status=gerr)
		if gerr ne 1 then begin
			print,'READG FAILED'
			equil={err:1}
			goto,GET_OUT
		endif
	endelse

    rhogrid=rho_rz(g,grid.r_grid/100.,grid.wc/100.,/do_linear)

	calculate_bfield,bp,br,bphi,bz,g
;	help,bp,br,bphi,bz
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
	ez=-(shift(epot,0,-1) - shift(epot,0,1))/(g.z[2]-g.z[0])
	
	;; Interpolate cylindrical fields onto (x,y,z) mesh then rotate vectors
	b=fltarr(3,grid.ng) & e=fltarr(3,grid.ng)
	for l=0l,grid.ng-1 do begin
		rgrid=(.01*grid.r_grid[l] - g.r[0])/(g.r[1]-g.r[0]) ; in grid units
		zgrid=(.01*grid.wc[l] - g.z[0])/(g.z[1]-g.z[0])    ; WWH 3/31/07
		bcylr=interpolate(br,[rgrid],[zgrid])
		ecylr=interpolate(er,[rgrid],[zgrid])
		bcylphi=interpolate(bphi,[rgrid],[zgrid])
		e[2,l]=interpolate(ez,[rgrid],[zgrid])
		b[2,l]=interpolate(bz,[rgrid],[zgrid])
		cph=cos(grid.phi_grid[l]) & sph=sin(grid.phi_grid[l])
		b[0,l]=(cph*bcylr - sph*bcylphi)
		b[1,l]=(sph*bcylr + cph*bcylphi)
		e[0,l]=cph*ecylr
		e[1,l]=sph*ecylr
	endfor

	if inputs.f90brems eq 0 then begin
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

	equil={g:g,rho_grid:rhogrid,rho_chords:rho_chords,b:b,e:e,err:0}
	GET_OUT:
	return,equil
END
