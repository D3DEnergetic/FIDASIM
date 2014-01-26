PRO read_inputs,file,inputs,save=save

    idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''
	
	if file_test(file) then begin
		openr,55,file
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
		readf,55,idum & load_fbm     = idum
		readf,55,idum & f90brems=idum
		readf,55,idum & calc_fida_wght    = idum
		readf,55,idum & calc_npa_wght    = idum
		readf,55,sdum;'# weight function settings:'
		readf,55,idum & ne_wght      = idum
		readf,55,idum & np_wght      = idum
		readf,55,idum & nphi_wght      = idum
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
        install_dir=+getenv('FIDASIM_DIR')
		inputs={instal_dir:install_dir,shot: shot, time: time,diag:diag,ps:ps $
	          , fidasim_runid:fidasim_runid $
	          , calc_fida_wght:calc_fida_wght,calc_npa_wght:calc_npa_wght,ne_wght:ne_wght,np_wght:np_wght,nphi_wght:nphi_wght $
              , ichan_wght:ichan_wght $
	          , nr_fast:nr_fast,nr_ndmc:nr_ndmc,nr_halo:nr_halo $
	          , calc_spec:calc_spec,npa:npa $
	          , einj:einj $
	          , load_neutrals:load_neutrals,load_fbm:load_fbm  $
	          , nx:nx, ny:ny, nz:nz $
	          , x0:xx[0],x1:xx[nx-1]+dx $
	          , y0:yy[0],y1:yy[ny-1]+dy $
	          , z0:zz[0],z1:zz[nz-1]+dz, isource:isource,origin:origin,alpha:alpha,beta:beta,err:0}
		if keyword_set(save) then save,inputs,filename='inputs.sav'
	endif else begin	
		inputs={err:1}
	endelse

END
PRO read_los,file,los,save=save

	if file_test(file) then begin
       	ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'Nchan',nchan
       	ncdf_varget,ncid,'xlens',xlens
       	ncdf_varget,ncid,'ylens',ylens
       	ncdf_varget,ncid,'zlens',zlens
       	ncdf_varget,ncid,'xlos',xlos
       	ncdf_varget,ncid,'ylos',ylos
       	ncdf_varget,ncid,'zlos',zlos
       	ncdf_varget,ncid,'ra',ra
       	ncdf_varget,ncid,'rd',rd
       	ncdf_varget,ncid,'sigma_pi',sigma_pi
       	ncdf_varget,ncid,'h',h
       	ncdf_varget,ncid,'los_wght',weight
       	ncdf_varget,ncid,'chan_id',chan_id

       	ncdf_close,ncid
		xyzlens=[[xlens],[ylens],[zlens]]
		xyzlos=[[xlos],[ylos],[zlos]]
		rlos=sqrt(xyzlos[*,0]^2+xyzlos[*,1]^2)
		
		los={nchan:nchan,rlos:rlos,xyzlens:xyzlens,xyzlos:xyzlos,$
			 ra:ra,rd:rd,h:h,chan_id:chan_id,sigma_pi:sigma_pi,weight:weight,err:0}
		if keyword_set(save) then save,los,filename='los.sav'
	endif else begin
		los={err:1}
	endelse
END
PRO read_fbm,file,FBM_struct,save=save,pitch_sign_convention=pitch_sign_convention
  ;---------------------------------------------------
  ; LOAD TRANSP fast ion distribution function
  ;---------------------------------------------------

  if not keyword_set(pitch_sign_convention) then pitch_sign_convention=-1.0

  ;;READ IN FILE
  if FILE_TEST(file) eq 1 then begin
	ncid=ncdf_open(file,/nowrite)

	ncdf_varget,ncid,'FBM_time',time
	ncdf_varget,ncid,'FBM_Nenergy',nenergy
	ncdf_varget,ncid,'FBM_Npitch',npitch
	ncdf_varget,ncid,'FBM_Ngrid',nzones
	ncdf_varget,ncid,'FBM_r2d',r2d
	ncdf_varget,ncid,'FBM_z2d',z2d
	ncdf_varget,ncid,'FBM_bmvol',bmvol
	ncdf_varget,ncid,'FBM',FBM
	ncdf_varget,ncid,'FBM_emin',emin
	ncdf_varget,ncid,'FBM_emax',emax
	ncdf_varget,ncid,'FBM_energy',energy
	ncdf_varget,ncid,'FBM_pmin',pmin
	ncdf_varget,ncid,'FBM_pmax',pmax
	ncdf_varget,ncid,'FBM_pitch',pitch
	ncdf_close,ncid

    ;in fidasim cdf file, pitch is defined in comp. to B-field.
    ;now revert to more common convention (to current) -> - sign.
    pitch *= pitch_sign_convention
    FBM_struct = { nenergy:nenergy, npitch:npitch, nzones:nzones $
                    ,FBM:FBM,energy:energy,pitch:pitch            $
                    ,emin:emin,emax:emax,pmin:pmin,pmax:pmax      $
                    ,r2d:r2d,z2d:z2d,bmvol:bmvol $
                    ,time:time $
                    ,pitch_sign_convention:pitch_sign_convention,err:0}
	if keyword_set(save) then begin
		fbm=FBM_struct
		save,fbm,filename='fbm.sav'
	endif
  endif else begin
	FBM_struct={err:1}
  endelse

END

PRO read_birth,file,birth,save=save

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'birth_dens',birth_dens
		ncdf_close,ncid
	
		birth={shot:shot,time:time,birth_dens:birth_dens,err:0}
		
		if keyword_set(save) then save,birth,filename='birth.sav'
	endif else begin
		birth={err:1}
	endelse

END

PRO read_grid,file,grid,save=save

	if file_test(file) then begin
	    ncid=ncdf_open(file,/nowrite)
	    ncdf_varget,ncid,'Nx',nx
	    ncdf_varget,ncid,'Ny',ny
	    ncdf_varget,ncid,'Nz',nz
		ncdf_varget,ncid,'u_grid',u_grid
		ncdf_varget,ncid,'v_grid',v_grid
		ncdf_varget,ncid,'w_grid',w_grid
		ncdf_varget,ncid,'x_grid',x_grid
		ncdf_varget,ncid,'y_grid',y_grid
		ncdf_varget,ncid,'z_grid',z_grid
		ncdf_varget,ncid,'r_grid',r_grid
		ncdf_varget,ncid,'phi_grid',phi_grid

	    ncdf_close,ncid

		grid={nx:nx,ny:ny,nz:nz,u_grid:u_grid,v_grid:v_grid,w_grid:w_grid,$
			  r_grid:r_grid,phi_grid:phi_grid,x_grid:x_grid,y_grid:y_grid,z_grid:z_grid,err:0}
		if keyword_set(save) then save,grid,filename='grid.sav'
	endif else begin
		grid={err:1}
	endelse

END
PRO read_spectra,nbi_halo_file,fida_file,spectra,save=save

    if file_test(nbi_halo_file) or file_test(fida_file) then begin
        nbi_halo_err=0
        fida_err=0
	    if file_test(nbi_halo_file) then begin
		    ncid=ncdf_open(nbi_halo_file,/nowrite)
		    ncdf_varget,ncid,'shot',shot
		    ncdf_varget,ncid,'time',time
		    ncdf_varget,ncid,'lambda',lambda
		    ncdf_varget,ncid,'full',full
		    ncdf_varget,ncid,'half',half
		    ncdf_varget,ncid,'third',third
		    ncdf_varget,ncid,'halo',halo
		    ncdf_varget,ncid,'brems',brems
		    ncdf_close,ncid
	    endif else begin
            full=0 & half=0 & third=0 & halo=0 & brems=0
	    endelse

	    if file_test(fida_file) then begin
       	    ncid=ncdf_open(fida_file,/nowrite)
       	    ncdf_varget,ncid,'shot',shot
       	    ncdf_varget,ncid,'time',time
       	    ncdf_varget,ncid,'lambda',lambda
       	    ncdf_varget,ncid,'spectra',fida
            ncdf_close,ncid
	    endif else begin
            fida=0
    	endelse
        spec=full+half+third+halo+fida+brems
	    spectra={shot:shot,time:time,lambda:lambda,spectra:spec,full:full,half:half,third:third,halo:halo,fida:fida,brems:brems,err:0}
    endif else spectra={err:1}
END

PRO read_neutrals,file,neutrals,save=save

	if file_test(file) then begin
		ncid=NCDF_OPEN(file,/nowrite)
		
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'fdens',fdens
		ncdf_varget,ncid,'hdens',hdens
		ncdf_varget,ncid,'tdens',tdens
		ncdf_varget,ncid,'halodens',halodens
		NCDF_CLOSE,ncid

     	neutrals={shot:shot,time:time,fdens:fdens, hdens:hdens,tdens:tdens, halodens:halodens,err:0}
		if keyword_set(save) then save,neutrals,filename='neutrals.sav'
	endif else begin
		neutrals={err:1}
  	endelse

END
PRO read_npa,file,npa,save=save

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'ipos',npaipos
		ncdf_varget,ncid,'fpos',npafpos
		ncdf_varget,ncid,'v',npav
		ncdf_varget,ncid,'wght',npawght
		ncdf_varget,ncid,'flux',flux
		ncdf_varget,ncid,'energy',energy
		ncdf_close,ncid

		if n_elements(npawght) gt 0 then begin
			npa={npaipos:npaipos,npafpos:npafpos,npav:npav,npawght:npawght,energy:energy,flux:flux,err:0}
			if keyword_set(save) then save,npa,filename='npa.sav'
		endif else begin
			print, 'attention, no particle reached the detector!'
			npa={err:1}
		endelse
	endif else begin
		npa={err:1}
	endelse
END
PRO read_plasma,plasma_file,plasma,save=save

	if file_test(plasma_file) then begin

       	ncid=ncdf_open(plasma_file,/nowrite)
       	ncdf_varget,ncid,'ti',ti
       	ncdf_varget,ncid,'te',te
       	ncdf_varget,ncid,'dene',dene
       	ncdf_varget,ncid,'deni',deni
       	ncdf_varget,ncid,'denp',denp
       	ncdf_varget,ncid,'denf',denf
       	ncdf_varget,ncid,'vrotx',vrotx
       	ncdf_varget,ncid,'vroty',vroty
       	ncdf_varget,ncid,'vrotz',vrotz
       	ncdf_varget,ncid,'zeff',zeff
       	ncdf_varget,ncid,'bx',bx
       	ncdf_varget,ncid,'by',by
       	ncdf_varget,ncid,'bz',bz
       	ncdf_varget,ncid,'ex',ex
       	ncdf_varget,ncid,'ey',ey
       	ncdf_varget,ncid,'ez',ez
       	ncdf_varget,ncid,'rho_grid',rho_grid

       	ncdf_close,ncid	
		
		plasma={rho_grid:rho_grid,te:te,ti:ti,$
				dene:dene,deni:deni,denp:denp,denf:denf,zeff:zeff,$
				vrotx:vrotx,vroty:vroty,vrotz:vrotz,bx:bx,by:by,bz:bz,$
				ex:ex,ey:ey,ez:ez,err:0}
		if keyword_set(save) then save,plasma,filename='plasma.sav'
	endif else begin
		plasma={err:1}
	endelse

END

PRO read_fida_weights,file,weights,save=save,pitch_sign_convention=pitch_sign_convention
	
	if not keyword_set(pitch_sign_convention) then pitch_sign_convention=-1.0

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'lambda',central_wavel
		ncdf_varget,ncid,'energy',energyarr
		ncdf_varget,ncid,'pitch',pitcharr
		ncdf_varget,ncid,'radius',radarr
		ncdf_varget,ncid,'theta',theta_arr
		ncdf_varget,ncid,'wfunct',weight_tot
		ncdf_close,ncid
		
		dE=abs(energyarr[1]-energyarr[0])
		dpitch=abs(pitcharr[1]-pitcharr[0])
		dwav=abs(central_wavel[1]-central_wavel[0])

		nwav=n_elements(central_wavel)
	    nen=n_elements(energyarr)
		npitch=n_elements(pitcharr)

	    pitcharr*=pitch_sign_convention
    	emax=max(energyarr)+0.5*dE

    	weights={shot:shot,time:time,nchan:n_elements(radarr),nen:nen $
            ,dE:dE,emax:emax,emin:0.,npitch:npitch,dpitch:dpitch   $
            ,nwav:nwav,dwav:dwav   $
            ,central_wavel:central_wavel,energyarr:energyarr,pitcharr:pitcharr $
            ,weight_tot:weight_tot   $
            ,angle:theta_arr,radius:radarr,err:0}
		if keyword_set(save) then save,weights,filename='fida_weights.sav'
	endif else begin
		weights={err:1}
	endelse

END

PRO read_npa_weights,file,weights,save=save,pitch_sign_convention=pitch_sign_convention
	
	if not keyword_set(pitch_sign_convention) then pitch_sign_convention=-1.0

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'energy',energyarr
		ncdf_varget,ncid,'pitch',pitcharr
		ncdf_varget,ncid,'radius',radarr
		ncdf_varget,ncid,'wfunct',weight_tot
        ncdf_varget,ncid,'flux',flux
		ncdf_close,ncid
		
		dE=abs(energyarr[1]-energyarr[0])
		dpitch=abs(pitcharr[1]-pitcharr[0])

	    nen=n_elements(energyarr)
		npitch=n_elements(pitcharr)

	    pitcharr*=pitch_sign_convention
    	emax=max(energyarr)+0.5*dE

    	weights={shot:shot,time:time,nchan:n_elements(radarr),nen:nen $
            ,dE:dE,emax:emax,emin:0.,npitch:npitch,dpitch:dpitch   $
            ,energyarr:energyarr,pitcharr:pitcharr $
            ,weight_tot:weight_tot,flux:flux   $
            ,radius:radarr,err:0}
		if keyword_set(save) then save,weights,filename='npa_weights.sav'
	endif else begin
		weights={err:1}
	endelse

END

PRO load_results,result_dir,results,save=save

	;;CHECK FOR SLASH
	slash=strmid(result_dir,0,1,/reverse_offset)
	if slash ne '/' then result_dir+='/'

	;;READ INPUTS
	read_inputs,result_dir+'inputs.dat',inputs

	runid=inputs.fidasim_runid
	input_file=inputs.fidasim_runid+'_inputs.cdf'

	;;READ GRID
	read_grid,result_dir+input_file,grid

	;;READ LOS
	read_los,result_dir+input_file,los

	;;READ PLASMA 
	read_plasma,result_dir+input_file,plasma

	;;READ SPECTRA
	read_spectra,result_dir+runid+'_nbi_halo_spectra.cdf', result_dir+runid+'_fida_spectra.cdf',spectra

	;;READ NEUTRALS
	read_neutrals,result_dir+runid+'_neutrals.cdf',neutrals

	;;READ NPA
	read_npa,result_dir+runid+'_npa.cdf',npa

	;;READ FBM
	read_fbm,result_dir+input_file,fbm

	;;READ FIDA WEIGHT FUNCTIONS
	read_fida_weights,result_dir+runid+'_fida_weight_function.cdf',fida_weights

	;;READ NPA WEIGHT FUNCTIONS
	read_npa_weights,result_dir+runid+'_npa_weight_function.cdf',npa_weights

	;;READ BIRTH PROFILE
	read_birth,result_dir+runid+'_birth.cdf',birth

	results={inputs:inputs,grid:grid,los:los,plasma:plasma,$
			spectra:spectra,neutrals:neutrals,$
			npa:npa,fbm:fbm,fida_weights:fida_weights,npa_weights:npa_weights,birth:birth}
	if keyword_set(save) then begin
		save,inputs,grid,los,plasma,fida,nbi_halo,neutrals,$
		npa,fbm,weights,birth,filename=inputs.fidasim_runid+'_results.sav',/compress
	endif
END
