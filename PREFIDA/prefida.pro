PRO rotate_uvw,uvw,Arot,Brot,Crot,updown,xyz 	
	;;rotate uvz by alpa alround z axis
  	if updown lt 0 then qrz=MATRIX_MULTIPLY(Arot,uvw)
  	if updown ge 0 then qrz=MATRIX_MULTIPLY(Brot,uvw)
  	;; rotate qrz_ray by phi_box on xyz_coordinates
  	xyz=MATRIX_MULTIPLY(Crot,qrz)
END

PRO rotate_points,x,y,z,Arot,Brot,Crot,xp,yp,zp
	xyz=transpose([[x],[y],[z]])
    rotate_uvw,xyz,Arot,Brot,Crot,1,xyz_p
    xp=transpose(xyz_p[0,*])
    yp=transpose(xyz_p[1,*])
    zp=transpose(xyz_p[2,*])
END

PRO make_rot_mat,ALPHA,BETA,Arot,Brot,Crot
    zero=0.d0
    one=1.d0
    ;;transformation matrix to rotate on NBI box axis by ALPHA
    Arot=dblarr(3,3)
    Arot[0,0]= cos(BETA)   & Arot[0,1]= zero    & Arot[0,2]= sin(BETA)
    Arot[1,0]= zero        & Arot[1,1]= one     & Arot[1,2]= zero
    Arot[2,0]=-sin(BETA)   & Arot[2,1]= zero    & Arot[2,2]= cos(BETA)
    ;; transformation matrix to rotate in vertical direction by beta
    Brot=dblarr(3,3)
    Brot[0,0]= cos(BETA)   & Brot[0,1]= zero    & Brot[0,2]= sin(BETA)
    Brot[1,0]= zero        & Brot[1,1]= one     & Brot[1,2]= zero
    Brot[2,0]=-sin(BETA)   & Brot[2,1]= zero    & Brot[2,2]= cos(BETA)
	;; Arot and Brot are exactly the same. I dont know why Ben Gieger had it
	;; but it will stay for now
    ;; transformation matrix to rotate towards the xy-coordinate system
    Crot=dblarr(3,3)
    Crot[0,0]= cos(ALPHA) & Crot[0,1]=-sin(ALPHA) & Crot[0,2]= zero
    Crot[1,0]= sin(ALPHA) & Crot[1,1]= cos(ALPHA) & Crot[1,2]= zero
    Crot[2,0]= zero 	  & Crot[2,1]= zero	      & Crot[2,2]= one
END

PRO make_fida_grid,inputs,grid,err

	err=1

	if inputs.alpha gt 2*!DPI or inputs.beta gt 2*!DPI then begin
		print,'Angles must be in radians'
		goto, GET_OUT
	endif
	make_rot_mat,inputs.alpha,inputs.beta,Arot,Brot,Crot

	nx=inputs.nx
	ny=inputs.ny
	nz=inputs.nz
	xdim1=inputs.xdim1
	xdim2=inputs.xdim2
	ydim1=inputs.ydim1
	ydim2=inputs.ydim2
	zdim1=inputs.zdim1
	zdim2=inputs.zdim2

	;Basic grid points
  	dx= (xdim2-xdim1) / double(nx)
  	dy= (ydim2-ydim1) / double(ny)
  	dz= (zdim2-zdim1) / double(nz)
	ng=long(nx)*long(ny)*long(nz)     ;; nr of cells

  	;;cell borders
  	xx=dx*dindgen(nx)+xdim1
  	yy=dy*dindgen(ny)+ydim1
  	zz=dz*dindgen(nz)+zdim1

	dr=[dx,dy,dz]
	drmin=min(dr)
	dv=dx*dy*dz

	;;cell centers
	xxc=xx+0.5d0*dx & yyc=yy+0.5d0*dy & zzc=zz+0.5d0*dz
	;; Put the basic grid into 1D array (useful for libkk routines)
	x=dblarr(ng) & y=dblarr(ng) & z=dblarr(ng)

	for i=0L,nx-1 do for j=0L,ny-1 do for k=0L,nz-1 do begin
		l=i+nx*j+nx*ny*k
		x[l]=xx[i] & y[l]=yy[j] & z[l]=zz[k]
	end

	;; Make the corresponding grid center arrays
	xc=x+0.5d0*dx & yc=y+0.5d0*dy & zc=z+0.5d0*dz
	
	;;Rotate all grid points to machine coordinates
	rotate_points,x,y,z,Arot,Brot,Crot,u,v,w
	rotate_points,xc,yc,zc,Arot,Brot,Crot,uc,vc,wc
	
	;;Change origin for rotated points
	u+=inputs.origin[0] & uc+=inputs.origin[0]
	v+=inputs.origin[1] & vc+=inputs.origin[1] 
	w+=inputs.origin[2] & wc+=inputs.origin[2]
	
	r_grid=sqrt(uc^2+vc^2)
	phi_grid=atan(vc,uc)
;	plot,xc,yc,psym=2,color=0,background=255
;    uvw_ray=[-1.d0,0.d0,0.d0]*1000.
;    rotate_uvw,uvw_ray,Arot,Brot,Crot,1,xyz_ray
;    oplot,[inputs.origin[0],inputs.origin[0]+xyz_ray[0]] $
;         ,[inputs.origin[1],inputs.origin[1]+xyz_ray[1]],thick=2,color=230

	grid={nx:nx,ny:ny,nz:nz,x:x,y:y,z:z,xx:xx,yy:yy,zz:zz,xc:xc,yc:yc,zc:zc,xxc:xxc,yyc:yyc,zzc:zzc,$
			dx:dx,dy:dy,dz:dz,dr:dr,drmin:drmin,dv:dv,ng:ng,u:u,v:v,w:w,$
			uc:uc,vc:vc,wc:wc,r_grid:r_grid,phi_grid:phi_grid}
	err=0	
	GET_OUT:
END
	
PRO prepare_beam,inputs,nbi,nbgeom
	
	rot_mat={err:1}
	isource=inputs.isource

	if nbi.pinj le 0. then begin
		print, 'the selected source nr',isource,' is not on!'
		goto, GET_OUT
	endif
	uvw_src=nbi.xyz_src-inputs.origin
	uvw_pos=nbi.xyz_pos-inputs.origin

	make_rot_mat,-inputs.alpha,inputs.beta,Drot,Erot,Frot
    rotate_uvw,uvw_src,Drot,Erot,Frot,1,xyz_src ;;rotate from machine to beam coordinates
    rotate_uvw,uvw_pos,Drot,Erot,Frot,1,xyz_pos
	
	xs=xyz_src[0] & ys=xyz_src[1] & zs=xyz_src[2]
	xp=xyz_pos[0] & yp=xyz_pos[1] & zp=xyz_pos[2]

	dis=sqrt( (xs-xp)^2.0d +(ys-yp)^2.0d + (zs-zp)^2.0d)
	BETA=double(asin((zp-zs)/dis))
	ALPHA=double(atan((yp-ys),(xp-xs))-!DPI)
	print,'BEAM CROSSOVER POINT:'
	print,xyz_pos
	print,'BEAM ROTATION ANGLES AS DEFINED BY fidasim.f90'
	print,'ALPHA: '
	print,ALPHA,FORMAT='(F20.10)'
	print,'BETA:'
	print,BETA,FORMAT='(F20.10)'

	;;MAKE ROTATION MATRICES 
	make_rot_mat,ALPHA,BETA,Arot,Brot,Crot

	nbgeom={isource:isource,alpha:ALPHA,beta:BETA,Arot:Arot,Brot:Brot,Crot:Crot,xyz_src:xyz_src,xyz_pos:xyz_pos,err:0}
	GET_OUT:
END

PRO los_track,coords,xyz_los_vec,xyspt,tcell,cell,ncell
	ri=xyspt
	vn=xyz_los_vec
	;; zeros are not good!
	index=where(vn eq 0,nind)
	if nind gt 0 then vn[index]=0.0001
	p=[0,0,0]
	dummy = min( abs( coords.xxc - ri[0] ), index )
	p[0]=index
	dummy = min( abs( coords.yyc - ri[1] ), index ) 
	p[1]=index
	dummy = min( abs( coords.zzc - ri[2] ), index ) 
	p[2]=index

	tcellh=fltarr(1000)
	cellh=fltarr(3,1000)
	m=0
	cellh[*,m]=p
	;;loop along line of sight
	while(m lt 1000) do begin
		l=p
		index=where(vn gt 0.d0)
		if index[0] ne -1 then begin
			l(index)=p(index)+1
		endif 
		if l[0] gt coords.nx-1 or l[1] gt coords.ny-1 or $
			l[2] gt coords.nz-1 then goto, out
		;;time needed to go into next cell
		dt_arr=fltarr(3)
		dt_arr[0] = ( coords.xx(l[0]) - ri[0] ) /vn[0]
		dt_arr[1] = ( coords.yy(l[1]) - ri[1] ) /vn[1]
		dt_arr[2] = ( coords.zz(l[2]) - ri[2] ) /vn[2] 
		dt=min(dt_arr,index)
 		ri[*] = ri[*] + vn[*]*dt 
		if vn[index] gt 0.d0 then begin	
			p[index]=p[index]+1
 		endif else begin
			p[index]=p[index]-1
		endelse
		if p[0] lt 0 or p[1] lt 0 or p[2] lt 0 then goto, out
		tcellh[m]=dt
		m=m+1
		cellh[*,m]=p[*] 
	endwhile
	out:
	;; Store results into compressed arrays!
	if m gt 1 then begin
		tcell= tcellh[0:m-2]
		cell = cellh[*,0:m-2]
	endif else begin
		tcell = -1
		cell = -1
	endelse
	ncell=m-1
END

PRO prepare_fida,inputs,grid,fida,chords

	nx=grid.nx
	ny=grid.ny
	nz=grid.nz
    ;;CALCULATE WEIGHTS
	err_arr=dblarr(fida.nchan)
	weight  = replicate(0.d0,nx,ny,nz,fida.nchan)
;	print, 'nchan:', fida.nchan

	make_rot_mat,-inputs.alpha,inputs.beta,Arot,Brot,Crot
	ulens=fida.xlens-inputs.origin[0] & ulos=fida.xlos-inputs.origin[0]
	vlens=fida.ylens-inputs.origin[1] & vlos=fida.ylos-inputs.origin[1]
	wlens=fida.zlens-inputs.origin[2] & wlos=fida.zlos-inputs.origin[2]

	rotate_points,ulens,vlens,wlens,Arot,Brot,Crot,xlens,ylens,zlens
	rotate_points,ulos,vlos,wlos,Arot,Brot,Crot,xlos,ylos,zlos

	for chan=0L, fida.nchan-1 do  begin
		xyzlens = [xlens[chan],ylens[chan],zlens[chan]]
        xyzlos  = [xlos[chan], ylos[chan], zlos[chan]]
;		print,xyzlos	
;		plot,[xlos[chan],xlens[chan]],[zlos[chan],zlens[chan]],color=0,background=255,xrange=[-600,-300],yrange=[-55.,55.]
;		oplot,grid.xc,grid.zc,psym=3,color=0
        vi    = xyzlos-xyzlens
        dummy = max(abs(vi),ic)
        nstep = fix(700./grid.dr[ic])
        vi    = vi/sqrt(vi[0]^2+vi[1]^2+vi[2]^2) ;; unit vector
;        if chan eq fida.nchan-1 then begin
;       	print, vi
;        endif
        xyz_pos = xyzlens
      ; find first grid cell
        for i=0L,nstep do begin
        	xyz_pos[0] = xyz_pos[0] + grid.dr[ic] * vi[0]/abs(vi[ic])
           	xyz_pos[1] = xyz_pos[1] + grid.dr[ic] * vi[1]/abs(vi[ic])
           	xyz_pos[2] = xyz_pos[2] + grid.dr[ic] * vi[2]/abs(vi[ic])
           	if xyz_pos[0] gt grid.xx[0] and xyz_pos[0] lt grid.xx[nx-1]+grid.dx and $ 
            	xyz_pos[1] gt grid.yy[0] and xyz_pos[1] lt grid.yy[ny-1]+grid.dy and $
              	xyz_pos[2] gt grid.zz[0] and xyz_pos[2] lt grid.zz[nz-1]+grid.dz then begin
              	goto, out
           	endif  
        endfor
        out:
      ; determine cells along the LOS
        if i lt nstep then begin
        	los_track,grid,vi,xyz_pos,tcell,cell,ncell
           	if ncell gt 1 then begin
            	for jj=0L,ncell-1 do begin
                	if finite(tcell[jj]) eq 0 then stop
                	;;  tcell is the length of the track (cm) as v is 1cm/s
                 	weight[cell[0,jj],cell[1,jj],cell[2,jj],chan]=tcell[jj]
              	endfor
           	endif else begin
             	print, 'WARNING: CHANNEL #'+strtrim(string(chan),1)+' ONLY CROSSES ONE CELL!'
				err_arr[chan]=1
           	endelse
        endif else begin
;        	print, 'LOS does not cross the simulation grid!'
;          	print,'chan: ', chan
			err_arr[chan]=1
        endelse
	endfor
    index=where(finite(weight) eq 0,nind)
    if nind gt 0 then begin
    	print,'weight set to 0. as it was NAN or Infinite!'
        weight[index]=0.
    endif	
	los=where(err_arr eq 0,nw)
	if nw eq 0 then begin
		print,'NO LINES OF SIGHT CROSSED THE SIMULATION GRID'
		err=1
	endif else begin
		print,strtrim(string(nw),1)+' OUT OF '+strtrim(string(fida.nchan),1)+' CHORDS CROSSED THE SIMULATION GRID'
		weight=weight[*,*,*,los]
		err=0
	endelse
	chords={nchan:fida.nchan,xlens:xlens,ylens:ylens,zlens:zlens,sigma_pi_ratio:fida.sigma_pi_ratio,$
			xlos:xlos,ylos:ylos,zlos:zlos,headsize:fida.headsize,los:los,weight:weight,err:err}
END

PRO transp_fbeam,inputs,grid,denf,err

	!P.charsize=1.
	!P.background=255 & !P.color=0
	;doplot=1
	;print, 'reading fast ion distribution function from transp output'
	cdftest=findfile(inputs.cdf_file)
	;print, '======================='
	if cdftest[0] eq '' then begin
		print,'ERROR: '+inputs.cdf_file+' WAS NOT FOUND'
		err=1
		goto,GET_OUT
	endif
	cdfid=NCDF_Open(inputs.cdf_file,/nowrite)
	;; Retrieve signals
	;; --------------------------------------------------------
	ncdf_varget, cdfid,'TRANSP_RUNID', runid
	ncdf_varget, cdfid,'TIME' , cdf_time    
	ncdf_varget, cdfid,'R2D'  , r2d     ; rposition of cells
	ncdf_varget, cdfid,'Z2D'  , z2d     ;zposition of cells
	ncdf_varget, cdfid,'BMVOL', bmvol   ; plasma volume
	ncdf_varget, cdfid,'E_D_NBI', energy ; central values
	ncdf_varget, cdfid,'A_D_NBI', pitch  ; central values
	ncdf_varget, cdfid,'F_D_NBI', FBM    ; fast-ion distribution function
	ncdf_varget, cdfid,'NTOT_D_NBI',ntot ; total number of fast ions
	ncdf_varget, cdfid,'RSURF', rsurf    ; flux surface
	ncdf_varget, cdfid,'ZSURF', zsurf    ; flux surface
	NCDF_Close,cdfid
	;; ----------- Check the time
	if abs(inputs.time-cdf_time) gt 0.02 then begin
		print, ' CDF file time:',cdf_time 
		print, 'WARNING: Time of CDF file and simulation disagree!'
	endif     
	;;-------------Convert eV-> keV
	energy=energy*1.0d-3          ;; fidasim needs energy in kev  
	fbm=fbm*1.0d3                 ;; now, this needs to be corrected
	;; as we now calculate with fast-ions/omega/keV/cm^3 
	;;------------Convert d_omega --> pitch
	;; Fast-ion distribution is given as a function of cm^3, energy
	;; and d_omega/4Pi. omega is the solild angle in 3D velocity space. In
	;; order to transform this to a function depending on pitch instead
	;; of d_omega/4PI, one has to multiply by 0.5!
	fbm=fbm*0.5
	;;------------Cut off energy
	;print,'EMIN: ',inputs.emin
	;print,'EMAX: ',inputs.emax
	w=where(energy ge inputs.emin and energy le inputs.emax)
	energy=energy[w]     
	fbm=fbm[w,*,*]
	index=where(fbm lt 0.,nind) 
	if nind gt 0. then fbm[index]=0.d0
	
	;;-----------store distribution funciton in binary file
	file =inputs.result_dir+inputs.runid+'/transp_fbm.bin'
	sz=size(FBM)
	nenergy=sz(1)
	npitch=sz(2)
	ngrid=sz(3)
	openw, lun, file, /get_lun
	writeu,lun , long(ngrid)
	for i=0L,ngrid-1 do writeu ,lun , float(r2d[i])
	for i=0L,ngrid-1 do writeu ,lun , float(z2d[i])
	writeu,lun , long(nenergy)
	writeu,lun , double(inputs.emin)
	writeu,lun , double(inputs.emax)
	for i=0L,nenergy-1 do writeu ,lun , double(energy[i])
	writeu,lun , long(npitch)
	writeu,lun , double(pitch[0])
	writeu,lun , double(pitch[npitch-1])
	for i=0L,npitch-1 do writeu ,lun , double(pitch[i])
	for i=0L,nenergy-1 do begin
		for j=0L,npitch-1 do begin
			for k=0L,ngrid-1 do begin    
				writeu ,lun , float(FBM[i,j,k]/max(FBM[*,*,k])) ;; normalized
			endfor
		endfor
	endfor
	close,lun
	free_lun, lun
	print, 'TRANSP distribution stored in BINARY: '+file
	
	;;----------Determine fast-ion density averaged over pitch and energy
	dE      = energy[2] - energy[1]
	dpitch  = pitch[2]  - pitch[1]
	fdens=total(reform(total(fbm,1)),1)*dE*dpitch
	;help,fdens,r2d,z2d
	;plot,r2d,fdens,color=0,background=255 
	;; ------map fdens on FIDASIM grid and sort out
	;; ------points outside the separatrix
	;;------Determine FIDAsim grid to map distribution function
	r_grid=fltarr(inputs.nx,inputs.ny,inputs.nz)
	z_grid=r_grid
	for i=0L,inputs.nx-1 do begin
		for j=0L,inputs.ny-1 do begin
			for k=0L,inputs.nz-1 do begin
				jj=i+inputs.nx*j+inputs.nx*inputs.ny*k
				r_grid[i,j,k]=grid.r_grid[jj]
				;z_grid[i,j,k]=grid.wc[jj]
			endfor      
		endfor
	endfor 
	f3d=dblarr(inputs.nx,inputs.ny,inputs.nz)*0.d0
	if ngrid le 220 then width=6. else width=4.
	for j=0L,inputs.ny-1 do begin
		rout=reform(r_grid[*,j,0])
		zout=grid.zzc[*]
		;zout=reform(z_grid[i,j,*])
		TRIANGULATE, r2d, z2d, tr     
		fdens2=griddata(r2d,z2d,fdens,xout=rout,yout=zout,/grid,/SHEPARDS,triangles=tr)
		for i=0L,inputs.nx-1 do begin
			for k=0L,inputs.nz-1 do begin
				r=rout[i]
				z=zout[k]
				a=sqrt((z2d-z)^2+(r2d-r)^2)
				amin=min(a)
				;; only write fdens2 if it is close to r2d,z2d grid
				if amin le width then f3d[i,j,k]=fdens2[i,k]         
			endfor
		endfor 
	endfor
	denf=dblarr(grid.ng)*0.d0   
	for i=0L,inputs.nx-1 do begin
		for j=0L,inputs.ny-1 do begin
			for k=0L,inputs.nz-1 do begin
				jj=i+inputs.nx*j+inputs.nx*inputs.ny*k
				denf[jj]=f3d[i,j,k] >0.
			endfor      
		endfor
	endfor 
	err=0
	GET_OUT:
END

PRO map_profiles,inputs,grid,equil,profiles,plasma,err
	;;-------------------------------------------------
	;; MAP kinetic profiles on FIDASIM grid
	;;------------------------------------------------- 

	rhomax=max(profiles.rho)
	ww=where(equil.rho_grid gt rhomax,nww)
	;;Electron density
	dene = 1.d-6 * interpol(profiles.dene,profiles.rho,equil.rho_grid) > 0. ;[1/cm^3]
	dene[ww]=0.001*1d13

	;;Zeff
	zeff = interpol(profiles.zeff,profiles.rho,equil.rho_grid) > 1.0
;	zeff[ww]=profiles.zeff[-1]
	zeff[ww]=1.0

	;;Impurity density
	deni = (zeff-1.)/(inputs.impurity_charge*(inputs.impurity_charge-1))*dene
	
	;;Proton density
	denp = dene-inputs.impurity_charge*deni
	print,total(deni)/total(denp)*100. ,' percent of impurities'
	
	;;Fast-ion density
	if keyword_set(inputs.nofida) then begin
    	denf=dene*0.d0
  	endif else begin
     	transp_fbeam,inputs,grid,denf,terr
		if terr eq 1 then begin
			print,'ERROR: FAILED TO MAP FAST ION DENSITY'
			err=1
			goto,GET_OUT
		endif
  	endelse
 
	;;Electron temperature
	te = 1.d-3 * interpol(profiles.te,profiles.rho,equil.rho_grid) > 0.001 ;keV
	te[ww]=0.001
	
	;;Ion temperature   
	ti = 1.d-3 * interpol(profiles.ti,profiles.rho,equil.rho_grid) > 0.001 ;keV
	if max(ti) gt 10. or max(te) gt 10. then begin
		print, 'WARNING:'
		print, 'Electron or Ion temperature greater than 10 keV'
		print, 'Look at the tables, they might only consider'
		print, 'temperatures less than 10keV!'
	endif
	ti[ww]=0.001

	;;Plasma rotation	
	vtor      =   interpol(profiles.vtor,profiles.rho,equil.rho_grid)*grid.r_grid ; [cm/s]  
	vtor[ww]  =   replicate(0.0,nww)*grid.r_grid[ww]
	vrot      =   fltarr(3,grid.ng)
	vrot[0,*] = - sin(grid.phi_grid)*vtor 
	vrot[1,*] =   cos(grid.phi_grid)*vtor
	vrot[2,*] =   0.d0 

	;;Rotate vector quantities to beam coordinates 
	make_rot_mat,-inputs.alpha,inputs.beta,Arot,Brot,Crot
	rotate_uvw,vrot,Arot,Brot,Crot,1,vrot_xyz ;;machine basis to beam basis
	rotate_uvw,equil.b,Arot,Brot,Crot,1,b_xyz
	rotate_uvw,equil.e,Arot,Brot,Crot,1,e_xyz

	;; test if there are NANs or Infinites in the input profiels
	index=where(finite([ti,te,dene,denp,zeff,denp,deni]) eq 0,nind)
	if nind gt 0 then stop
	;;-------SAVE-------
	plasma={rho_grid:equil.rho_grid, b:b_xyz,e:e_xyz,ab:inputs.ab,ai:inputs.ai,te:te, $
			ti:ti,vtor:vtor,vrot:vrot_xyz,vrot_uvw:vrot,dene:dene,denp:denp,deni:deni,denf:denf $
			,zeff:zeff}
	err=0
	GET_OUT: 
END

FUNCTION sinterpol,v,x,u,sortt=sortt,_extra=_extra
	if n_elements(sortt) lt 1 then sortt=0

	if sortt then begin
		ind=sort(X)
	endif else begin
		ind=lindgen(n_elements(x))
	endelse


	return,interpol(v[ind],x[ind],u,_extra=_extra)
END

PRO brems,result_dir,det,profiles,equil
	;; Calculates visible bremsstrahlung along FIDA sightlines
	;; WWH 6/2013

	;; INPUT
	;; result_dir directory to write output
	;; det	     structure with detector lines of sight
	;; profiles   structure with plasma profiles vs. rho

	;; OUTPUT
	;; file with the surface radiance (ph/s-m2-nm-sr) for
	;; each sightline

	;*******************************************************************
	;*****************************************
	; Plasma parameters
	rho=profiles.rho
	te=profiles.te              ; eV
	dene=profiles.dene*1.e-6    ; cm^-3
	zeff=profiles.zeff
	
	; Require non-zero values for te and ne
	w=where(te le 0. or dene le 0.,nw)
	if nw gt 0 then begin
		rho=rho[0:w[0]-1]
		te=te[0:w[0]-1]
		dene=dene[0:w[0]-1]
		zeff=zeff[0:w[0]-1]
	endif
	w=where(zeff lt 1.,nw) & if nw gt 0 then zeff[w]=1.
	w=where(zeff gt 6.,nw) & if nw gt 0 then zeff[w]=6.
	rhomax=max(rho,nr) & nr+=1

	;**********************************************
	; Constants in calculation
	lambda=6561.	; average wavelength (Angstroms)
	h_planck=4.135667e-15  ; [eV/s]
	c0=2.9979e8 ; [m/s]
	
	; Visible bremsstrahlung emissivity versus rho
	gaunt=5.542-(3.108-alog(te/1000.))*(0.6905-0.1323/zeff)
	emisrho=10.*7.57d-9*gaunt*dene^2*zeff/(lambda*sqrt(te)) $
               *exp(-h_planck*c0/(lambda*te))

	;***********************************************
	;NOW do line integration to get surface radiance
	;***********************************************
	nchan=det.nchan
	vbline=replicate(0.,nchan)

	for i=0,nchan-1 do begin
		rhospath=equil.rho_chords.rhos[*,i]
		rhospath=rhospath[where(finite(rhospath))]
        vbepath=sinterpol(emisrho,rho,rhospath,/sort)
        wgtr1=where(rhospath ge rhomax,nwgtr1)
        ;set emission at radii outside of max rho to zero
        if nwgtr1 ge 0 then vbepath[wgtr1]=0.0
		vbline[i]=total(vbepath)*equil.rho_chords.ds*0.01   ; (ph/s-m2-nm-sr)
	endfor  ; channel loop

	;--------------------------------
	; Save results

	print,vbline
	file =result_dir+'/bremsstrahlung.bin'
	openw, lun, file, /get_lun
	for i=0,nchan-1 do writeu,lun, double(vbline[i])
	close,lun
	free_lun,lun
	print, 'Bremsstrahlung stored in BINARY: '+file
END

PRO prefida,input_pro,plot=plot

	COMPILE_OPT DEFINT32

	if n_elements(input_pro) eq 0 then begin
		print,'NEEDS INPUT FILE'
		goto,GET_OUT
	endif

	;;CALL INPUT PROCEDURE/FILE
	CALL_PROCEDURE,input_pro,inputs

	;;MAKE DIRECTORIES IF THEY DONT EXIST
	if file_test(inputs.result_dir,/directory) eq 0 then begin
		spawn,'mkdir '+inputs.result_dir
	endif
	if file_test(inputs.result_dir+inputs.runid,/directory) eq 0 then begin
		spawn,'mkdir '+inputs.result_dir+inputs.runid
	endif

	;;ADD DEVICE DIRECTORY TO PATH
	!path = !path + ":" + expand_path("+"+inputs.install_dir+inputs.device)
	
	;;ADD INSTALL DIRECTORY TO PATH
	!path = !path + ":" + expand_path(inputs.install_dir)

	;;MAKE FIDA GRID
	make_fida_grid,inputs,grid,err
	if err eq 1 then begin
		print,'GRID CREATION FAILED. EXITING...'
		goto,GET_OUT
	endif else err=0

	;;CALL DEVICE ROUTINES THAT GET BEAM GEOMETRY, FIDA DIAGNOSTIC INFO, PROFILES, and the grid in flux coord.
	CALL_PROCEDURE, strlowcase(inputs.device)+'_routines',inputs,grid, nbi, fida, profiles, equil,err
	if err eq 1 then begin
		print, 'DEVICE ROUTINES FAILED. EXITING...'
		goto,GET_OUT
	endif else err=0

	;;BEAM PRE PROCESSING
	prepare_beam,inputs,nbi,nbgeom
	if nbgeom.err eq 1 then begin
		print,'BEAM PREPROCESSING FAILED. EXITING...'
		goto, GET_OUT
	endif else err=0

	;;FIDA PRE PROCESSING 
	if inputs.no_spectra ne 1 then begin
		prepare_fida,inputs,grid,fida,chords
		if chords.err eq 1 then begin
			print,'FIDA PREPROCESSING FAILED. EXITING...'
			goto, GET_OUT
		endif
	endif else begin
		err=0
		weight=replicate(0.d0,inputs.nx,inputs.ny,inputs.nz,1)
		chords={weight:weight,err:err}
	endelse

	;;MAP PROFILES ONTO GRID
	map_profiles,inputs,grid,equil,profiles,plasma,err
	if err eq 1 then begin
		print,'PROFILE MAPPING FAILED'
		goto,GET_OUT
	endif else err=0

    ;; Calculate bremsstrahlung if desired
	if inputs.f90brems eq 0 then $
		brems,inputs.result_dir+inputs.runid,fida,profiles,equil

    ;; Plot grid, beam, sightlines, and equilibrium
	if keyword_set(plot) then begin
		CALL_PROCEDURE, strlowcase(inputs.device)+'_plots',inputs,grid, nbi, fida, equil,nbgeom,plasma
	endif

	;;WRITE FIDASIM INPUT FILES
	file = inputs.result_dir+inputs.runid+'/inputs.dat'
	openw, 55, file
	printf,55,'# FIDASIM input file created: ', systime()
	printf,55, inputs.install_dir
	printf,55, inputs.shot         ,f='(i6,"         # shotnumber")'  
	printf,55, inputs.time,f='(1f8.5,"       # time")'
	printf,55, inputs.runid
	printf,55,' ',inputs.fida_diag, '           # diagnostic'
	printf,55,'# general settings:'
	printf,55, inputs.no_spectra,f='(i2,"             # no spectra")'
	printf,55, inputs.nofida,f='(i2,"             # only NBI+HALO")'
	printf,55, inputs.npa          ,f='(i2,"             # NPA simulation")'
	printf,55, inputs.load_neutrals,f='(i2,"             # load NBI+HALO density")'
    printf,55, inputs.guidingcenter,f='(i2,"             # 0 for full-orbit F")'
    printf,55, inputs.f90brems,f='(i2,"             # 0 reads IDL v.b.")'
	printf,55, inputs.calc_wght,f='(i2,"             # calculate wght function")'
	printf,55,'# weight function settings:'
	printf,55, inputs.nr_wght,f='(i9,"      # number velocities")'
	printf,55, inputs.ichan_wght,f='(i3,"      # channel for weight function")'
	printf,55, inputs.emax_wght,f='(1f12.2,"       # emax for weights")'
	printf,55, inputs.dwav_wght,f='(1f12.5,"       # dwav")'
	printf,55, inputs.wavel_start_wght,f='(1f12.5,"       # wavel_start")'
	printf,55, inputs.wavel_end_wght,f='(1f12.5,"       # wavel_end")'
	printf,55,'# Monte Carlo settings:'
	printf,55, inputs.nr_fida,f='(i9,"      # number of FIDA mc particles")'  
	printf,55, inputs.nr_ndmc,f='(i9,"      # number of NBI mc particles")' 
	printf,55, inputs.nr_halo,f='(i9,"      # number of HALO mc particles")'
	printf,55, inputs.impurity_charge,f='(i2,"             # Impurity charge")'
	printf,55,'# Location of transp cdf file:'
	printf,55, inputs.cdf_file
	printf,55,'# discharge parameters:'
	printf,55, inputs.btipsign,f='(i3,"            # B*Ip sign")'
	printf,55, inputs.ai,f='(1f7.4,"        # plasma mass")'
	printf,55, inputs.ab,f='(1f7.4,"        # NBI mass")'
	printf,55,'# wavelength grid:'
	printf,55, inputs.nlambda,f='(1i5,"          # nlambda")'
	printf,55, inputs.lambdamin,f='(1f9.3,"      # lambda min")'
	printf,55, inputs.lambdamax,f='(1f9.3,"      # lambda max")'
	printf,55,'# simulation grid: '
	printf,55, inputs.nx,f='(1i3,"            # nx")'
	printf,55, inputs.ny,f='(1i3,"            # ny")'
	printf,55, inputs.nz,f='(1i3,"            # nz")'  
	for i=0L,inputs.nx-1 do begin   ;; cell borders          
		printf,55,grid.xx[i],f='(1f9.4,"      # xx[i]")'
	endfor
	for i=0L,inputs.ny-1 do begin
		printf,55,grid.yy[i],f='(1f9.4,"      # yy[i]")'
	endfor
	for i=0L,inputs.nz-1 do begin
		printf,55,grid.zz[i],f='(1f9.4,"      # zz[i]")'
	endfor
	printf,55,'# Neutral beam injection:'
	printf,55, nbi.BMWIDRA,f='(1f9.4,"      # NBI half width horizontal")'
	printf,55, nbi.BMWIDZA,f='(1f9.4,"      # NBI half width vertical")'
	ii=inputs.isource[0]
	printf,55, ii,f='(1i2,"             # Nr of NBI")'
	printf,55,nbi.divy[0],f='(1f10.7,"     #divergence y of full comp")'
	printf,55,nbi.divy[1],f='(1f10.7,"     #divergence y of half comp")'
	printf,55,nbi.divy[2],f='(1f10.7,"     #divergence y of third comp")'
	printf,55,nbi.divz[0],f='(1f10.7,"     #divergence z of full comp")'
	printf,55,nbi.divz[1],f='(1f10.7,"     #divergence z of half comp")'
	printf,55,nbi.divz[2],f='(1f10.7,"     #divergence z of third comp")'
	printf,55,nbi.focy,f='(1f10.2,"      # focal length in y")' 
	printf,55,nbi.focz,f='(1f10.2,"      # focal length in z")' 
	printf,55,nbi.einj,f='(1f9.4,"      # injected energy [keV]")' 
	printf,55,nbi.pinj,f='(1f9.4,"      # injected power [MW]")'  
	printf,55,'# Species-mix (Particles):'
	printf,55 ,nbi.full,f='(1f9.6,"      # full energy")' 
	printf,55 ,nbi.half,f='(1f9.6,"      # half energy")'  
	printf,55 ,nbi.third,f='(1f9.6,"      # third energy")' 
	printf,55, '#position of NBI source in xyz coords:'
	printf,55,nbgeom.xyz_src[0],f='(1f9.4,"      # x [cm]")' 
	printf,55,nbgeom.xyz_src[1],f='(1f9.4,"      # y [cm]")' 
	printf,55,nbgeom.xyz_src[2],f='(1f9.4,"      # z [cm]")' 
	printf,55,'# 3 rotation matrizes 3x3'
	for j=0,2 do begin
		for k=0,2 do begin
			printf,55 ,nbgeom.Arot[j,k] ;; rotation in the top-down view plane
			printf,55 ,nbgeom.Brot[j,k] ;; vertical rotation
			printf,55 ,nbgeom.Crot[j,k] ;; vertical rotation
		endfor 
	endfor
	close,55
	print, 'Inputs stored in data file: '+file

	;;WRITE PLASMA PARAMETERS TO BINARY
	file =inputs.result_dir+inputs.runid+'/plasma.bin'
	openw, lun, file, /get_lun
	writeu,lun , long(inputs.nx)
	writeu,lun , long(inputs.ny)
	writeu,lun , long(inputs.nz)
	for ix=0L,inputs.nx-1 do begin
		for iy=0L,inputs.ny-1 do begin
			for iz=0L,inputs.nz-1 do begin
				i=ix+inputs.nx*iy+inputs.nx*inputs.ny*iz
				writeu,lun $
				, double(plasma.te[i])     , double(plasma.ti[i])    $
				, double(plasma.dene[i])   , double(plasma.denp[i])  $
				, double(plasma.deni[i])   , double(plasma.vrot[0,i])$
				, double(plasma.vrot[1,i]) , double(plasma.vrot[2,i])$
				, double(plasma.b[0,i])    , double(plasma.b[1,i])   $
				, double(plasma.b[2,i])    , double(plasma.e[0,i])   $
				, double(plasma.e[1,i])    , double(plasma.e[2,i])   $
				, double(plasma.rho_grid[i]),double(plasma.denf[i])  $
				, double(plasma.zeff[i])
			endfor
		endfor
	endfor
	close,lun
	free_lun, lun
	print, 'Plasma parameters stored in BINARY: '+file

	;;WRITE LINE OF SIGHT (LOS) INFORMATION TO BINARY
	if inputs.no_spectra ne 1 then begin
		file =inputs.result_dir+inputs.runid+'/los.bin'
		los=chords.los
		openw, lun, file, /get_lun
		writeu,lun , long(n_elements(los))
		for chan=0L,n_elements(los)-1 do begin
			writeu,lun, double(chords.xlens[los[chan]])
			writeu,lun, double(chords.ylens[los[chan]])
			writeu,lun, double(chords.zlens[los[chan]])
			writeu,lun, double(chords.headsize[los[chan]]) ;; headsize is used for NPA
			writeu,lun, double(chords.xlos[los[chan]])
			writeu,lun, double(chords.ylos[los[chan]])
			writeu,lun, double(chords.zlos[los[chan]])
		endfor
		writeu,lun , double(chords.sigma_pi_ratio)
		for i=0L,inputs.nx-1 do begin
			for j=0L,inputs.ny-1 do begin
				for k=0L,inputs.nz-1 do begin
					for chan=0L,n_elements(los)-1 do begin
						writeu ,lun , float(chords.weight[i,j,k,chan])
					endfor
				endfor
			endfor
		endfor
		close,lun
		free_lun, lun
		print, 'LOS parameters stored in BINARY: '+file
	endif

	print,''
	print,''
	print, 'To run FIDASIM use the following command'
	print, inputs.install_dir+'fidasim '+inputs.result_dir+inputs.runid
	print,''
	print,''
	GET_OUT:
END

