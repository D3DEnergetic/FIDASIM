PRO read_inputs,file,inputs

    idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''
	
	if file_test(file) then begin
		openr,55,file
		readf,55,sdum                 ;& print, sdum
		readf,55,sdum                 ;& print, sdum
		readf,55,idum & shot = long(idum) 
		readf,55,fdum & time = float(fdum)
		readf,55,sdum & fidasim_runid = sdum
		readf,55,sdum & diag          = strmid(sdum,1,3) 
		readf,55,sdum
		readf,55,idum & calc_spec    = idum
		readf,55,idum & ps           = idum
		readf,55,idum & npa          = idum      ;& print, 'npa: ', npa
		readf,55,idum & load_neutrals= idum
		readf,55,idum & f90brems=idum
		readf,55,idum & calc_wght    = idum
		readf,55,sdum;'# weight function settings:'
		readf,55,idum & nr_wght      = idum
		readf,55,idum & ichan_wght   = idum
		readf,55,sdum;inputs.emax_wght
		readf,55,sdum;inputs.dwav_wght
		readf,55,sdum;inputs.wavel_start_wght
		readf,55,sdum;inputs.wavel_end_wght
		readf,55,sdum;'# Monte Carlo settings:'
		readf,55,fdum & nr_fast = long(fdum)
		readf,55,fdum & nr_ndmc = long(fdum)
		readf,55,fdum & nr_halo = long(fdum)
		readf,55,sdum;impurity_charge,f='(i2,"             # Impurity charge")'
		readf,55,sdum;'# Location of transp cdf file:'
		readf,55,sdum & transp_runid=sdum
		readf,55,sdum                 ;& print, sdum
		readf,55,sdum                 ;& print, sdum
		readf,55,fdum & ai=fdum       ;& print, 'ai: ', ai
		readf,55,fdum & ab=fdum       ;& print, 'ab: ', ai
		for i=0,4 do readf,55,sdum
		origin=fltarr(3)
		for i=0,2 do begin 
			readf,55,fdum & origin[i]=fdum
		endfor
		readf,55,ddum & alpha=ddum
		readf,55,ddum & beta=ddum
		readf,55,idum & nx=idum       ;& print, 'nx: ', nx
		readf,55,idum & ny=idum       ;& print, 'ny: ', ny
		readf,55,idum & nz=idum       ;& print, 'nz: ', nz
		xx=fltarr(nx)
		yy=fltarr(ny)
		zz=fltarr(nz)
		for i=0,nx-1 do begin
			readf,55,fdum  & xx[i]=fdum	
		endfor  
		dx=xx[1]-xx[0]
		for i=0,ny-1 do begin
			readf,55, fdum  & yy[i]=fdum 
		endfor
		dy=yy[1]-yy[0]
		for i=0,nz-1 do begin
			readf,55,fdum   & zz[i]=fdum
		endfor
		dz=zz[1]-zz[0]

		;;READ NEUTRAL BEAM DATA
		readf,55,sdum
		readf,55,fdum & bmwidra=fdum
		readf,55,fdum & bmwidza=fdum
		nr_src=1
		xyz_src=fltarr(3)
		sadum = strarr(3)
		Arot   =fltarr(3,3)
		Brot   =fltarr(3,3)
		Crot   =fltarr(3,3)
		divy = fltarr(3,8)
		divz = fltarr(3,8)
		focy = fltarr(8)
		focz = fltarr(8)
		readf,55,fdum & isource=fdum
		readf,55,sadum & divy[*,isource] = double(sadum)
		readf,55,sadum & divz[*,isource] = double(sadum)
		readf,55,fdum & focy[isource]=fdum
		readf,55,fdum & focz[isource]=fdum
		readf,55,fdum & einj=fdum
		mass_u       = 1.6605402d-27                   ; [kg]
		e0           = 1.60217733d-19                  ; [C]             
		vinj= sqrt(2.d0*einj*1.e3*e0/(ab*mass_u))*1.d2 ;; [cm/s]       
		for i=0,5 do readf,55,sdum
			readf,55,fdum & xyz_src[0]=fdum 
			readf,55,fdum & xyz_src[1]=fdum
			readf,55,fdum & xyz_src[2]=fdum
			readf,55,sdum ;& print, sdum
		for j=0,2 do begin
			for k=0,2 do begin
		    	readf,55 ,fdum & Arot[j,k]=fdum
      			readf,55 ,fdum & Brot[j,k]=fdum
			    readf,55 ,fdum & Crot[j,k]=fdum
    	 	endfor 
		endfor
	
		close,55

		inputs={shot: shot, time: time,diag:diag,ps:ps $
	          , transp_runid:transp_runid $
	          , fidasim_runid:fidasim_runid $
	          , calc_wght:calc_wght,nr_wght:nr_wght,ichan_wght:ichan_wght $
	          , nr_fast:nr_fast,nr_ndmc:nr_ndmc,nr_halo:nr_halo $
	          , calc_spec:calc_spec,npa:npa $
	          , einj:einj $
	          , load_neutrals:load_neutrals  $
	          , nx:nx, ny:ny, nz:nz $
	          , x0:xx[0],x1:xx[nx-1]+dx $
	          , y0:yy[0],y1:yy[ny-1]+dy $
	          , z0:zz[0],z1:zz[nz-1]+dz, isource:isource,origin:origin,alpha:alpha,beta:beta,err:0}
	endif else begin	
		inputs={err:1}
	endelse

END
