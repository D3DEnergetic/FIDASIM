FUNCTION d3d_equil,inputs,grid

    g=readg(inputs.shot,inputs.time*1000,RUNID=inputs.equil,status=gerr)
    if gerr ne 1 then goto,GET_OUT
    rhogrid=rho_rz(g,grid.r_grid,grid.zc)

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
		zgrid=(.01*grid.zc[l] - g.z[0])/(g.z[1]-g.z[0])    ; WWH 3/31/07
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
	equil={rho_grid:rhogrid,b:b,e:e}
	GET_OUT:
	return,equil
END
