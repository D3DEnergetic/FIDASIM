;; --------------------------------------------------------
;; ------------------ LOAD_FIDASIM_RESULTS ----------------
;; --------------------------------------------------------
;; ROUTINE OF FIDASIM to read the results from binary and ascii files
;; written by the FORTRAN routine and by the inputs routine
;; (dalpha_inputs.pro). The output is a structure, containing all
;; information!
;; written by Benedikt Geiger and Markus Weiland 2013
;; additional subroutines
pro read_fbm,result_dir,fbm_struct,pitch_sign_convention
  ;---------------------------------------------------
  ; LOAD TRANSP fast ion distribution function
  ;---------------------------------------------------
  sdum=string(1,f='(i17)')
  ddum=1.d0
  idum=1L
  file =result_dir+'/transp_fbm.bin'
  FBM_struct = ptr_new()  ;If file is not availible, leave !null pointer.
  ;nur wenn File existiert:
  if FILE_TEST(file) eq 1 then begin
     openr, lun, file, /get_lun
     readu,lun, sdum & cdf_file=sdum
     readu,lun, ddum & time=ddum
     ;; SPATIAL GRID
     readu,lun , idum & nzones=idum
     r2d=dblarr(nzones)
     z2d=dblarr(nzones)
     bmvol=dblarr(nzones)
     readu, lun, r2d
     readu, lun, z2d
     readu, lun, bmvol
     ;; ENERGY GRID
     nenergy=1L
     emin   =1.d0 
     emax   =1.d0
     readu,lun , idum & nenergy=idum
     readu,lun , ddum & emin=ddum
     readu,lun , ddum & emax=ddum
     energy=dblarr(nenergy)
     readu, lun, energy
     ;; PITCH GRID
     readu,lun , idum & npitch=idum
     readu,lun,ddum & pmin=ddum
     readu,lun,ddum & pmax=ddum
     pitch=dblarr(npitch)
     readu,lun, pitch
     FBM = dblarr(nenergy,npitch,nzones)
     readu,lun, FBM
     close,lun
     free_lun, lun
     ;in fidasim binary file, pitch is defined in comp. to B-field.
     ;now revert to more common convention (to current) -> - sign.
     pitch *= pitch_sign_convention
     FBM_struct = { nenergy:nenergy, npitch:npitch, nzones:nzones $
                    ,FBM:FBM,energy:energy,pitch:pitch            $
                    ,emin:emin,emax:emax,pmin:pmin,pmax:pmax      $
                    ,r2d:r2d,z2d:z2d,bmvol:bmvol $
                    ,time:time,cdf_file:cdf_file $
                    ,pitch_sign_convention:pitch_sign_convention}
  endif
end

@LIB/set_det.pro
pro load_fidasim_results,results,result_dir_in $
                         ,only_inputs=only_inputs $
                         ,only_neutrals=only_neutrals, no_fbm=no_fbm $
                         ,only_fbm=only_fbm

  if n_elements(result_dir_in) eq 0 then begin
     result_dir_in = dialog_pickfile(/dir, path='~/FIDASIM/RESULTS')
  endif

  result_dir=result_dir_in
  print,'Loading fidasim results of run:'
  print, result_dir

  ;; define dummy variables
  string_dum=''
  idum=1L
  fdum=1.e0
  ddum=1.d0
  ;Pitch Sign Convention: 1 for FIDASIM-Convention (-1 == co-current)
  pitch_sign_convention = -1 


  ;this routine expects the path without the leading '/'.
  ;if its there - remove it:
  strlen = strlen(result_dir)
  if STRPOS(result_dir, '/', /REVERSE_SEARCH) eq strlen-1 then begin
     result_dir=strmid(result_dir, 0, strlen-1)
  endif
  

  if keyword_set(only_neutrals) then goto, only_neut

  ;; ------------------------------------------------------------
  ;; ----------------- LOAD THE INPUT FILE ----------------------
  ;; ------------------------------------------------------------
  input_file=result_dir+'/inputs.dat'  
  openr,55,input_file
  readf,55,string_dum                 ;& print, string_dum
  readf,55,string_dum                 ;& print, string_dum
  readf,55,idum & shot = long(idum) 
  readf,55,fdum & time = float(fdum)
  readf,55,string_dum & fidasim_runid = string_dum
  readf,55,string_dum & diag          = strmid(string_dum,1,3) 
  readf,55,string_dum
  readf,55,idum & calc_spec    = idum
  readf,55,idum & ps           = idum
  readf,55,idum & npa          = idum      ;& print, 'npa: ', npa
  readf,55,idum & load_neutrals= idum
  readf,55,idum & f90brems=idum
  readf,55,idum & calc_wght    = idum
  readf,55,string_dum;'# weight function settings:'
  readf,55,idum & nr_wght      = idum
  readf,55,idum & ichan_wght   = idum
  readf,55,string_dum;inputs.emax_wght
  readf,55,string_dum;inputs.dwav_wght
  readf,55,string_dum;inputs.wavel_start_wght
  readf,55,string_dum;inputs.wavel_end_wght
  readf,55,string_dum;'# Monte Carlo settings:'
  readf,55,fdum & nr_fast = long(fdum)
  readf,55,fdum & nr_ndmc = long(fdum)
  readf,55,fdum & nr_halo = long(fdum)
  readf,55,string_dum;impurity_charge,f='(i2,"             # Impurity charge")'
  readf,55,string_dum;'# Location of transp cdf file:'
  readf,55,string_dum & transp_runid=string_dum
  readf,55,string_dum                 ;& print, string_dum
  readf,55,string_dum                 ;& print, string_dum
  readf,55,fdum & ai=fdum       ;& print, 'ai: ', ai
  readf,55,fdum & ab=fdum       ;& print, 'ab: ', ai
  for i=0,4 do readf,55,string_dum
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
  ;read neutral beam data:
  readf,55,string_dum
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
  for i=0,5 do readf,55,string_dum
  readf,55,fdum & xyz_src[0]=fdum 
  readf,55,fdum & xyz_src[1]=fdum
  readf,55,fdum & xyz_src[2]=fdum
  readf,55,string_dum ;& print, string_dum
  for j=0,2 do begin
     for k=0,2 do begin
        readf,55 ,fdum & Arot[j,k]=fdum
        readf,55 ,fdum & Brot[j,k]=fdum
        readf,55 ,fdum & Crot[j,k]=fdum
     endfor 
  endfor
  if not eof(55) then begin    ;in earlier versions rhostring is not saved.
     readf,55,string_dum           
     readf,55,string_dum & rhostr=string_dum
     readf,55,string_dum           
     readf,55,fdum & rotate=fdum
  endif
  close,55
  if n_elements(rhostr) eq 0 then rhostr='rho'
            ;if rhostring is not saved, then just write rho ...
  if n_elements(rotate) eq 0 then rotate=0.
  equil_exp='AUGD'
  equil_diag='EQH'
  inputs={shot: shot, time: time,diag:diag,ps:ps $
          , equil_exp: equil_exp,equil_diag:equil_diag $
          , transp_runid:transp_runid $
          , fidasim_runid:fidasim_runid $
          , result_dir:result_dir $
          , calc_wght:calc_wght,nr_wght:nr_wght,ichan_wght:ichan_wght $         
          , nr_fast:nr_fast,nr_ndmc:nr_ndmc,nr_halo:nr_halo $
          , calc_spec:calc_spec,npa:npa $
          , einj:einj $
          , load_neutrals:load_neutrals  $
          , rotate:rotate  $
          , nx:nx, ny:ny, nz:nz $
          , x0:xx[0],x1:xx[nx-1]+dx $
          , y0:yy[0],y1:yy[ny-1]+dy $
          , z0:zz[0],z1:zz[nz-1]+dz, rhostr:rhostr, isource:isource,origin:origin,alpha:alpha,beta:beta}
  if keyword_set(only_inputs) then begin
     results={inputs:inputs}
     return
  endif


  ;; --------------------------------------------------------------------
  ;; -DEFINE ADDITIONAL structions with information from the input file -
  ;; --------------------------------------------------------------------
  ;nx=inputs.nx & ny=inputs.ny & nz=inputs.nz
  ng=long(nx)*long(ny)*long(nz)     ;; nr of cells
  ;dx=inputs.dx & dy=inputs.dy & dz=inputs.dz
  dr=[dx,dy,dz]        ;; size of cells
  drmin=min(dr)        ;; minimal size
  dv=dr[0]*dr[1]*dr[2] ;;volume
  xxc=xx+0.5d0*dx & yyc=yy+0.5d0*dy & zzc=zz+0.5d0*dz
  rrc_grid=dblarr(nx,ny,nz)
  phi_grid=dblarr(nx,ny,nz)
  zzc_grid=dblarr(nx,ny,nz)
  for i=0L,nx-1 do begin
     for j=0L,ny-1 do begin
        for k=0L,nz-1 do begin
           rrc_grid[i,j,k]=sqrt(xxc[i]^2+yyc[j]^2)
           phi_grid[i,j,k]=atan(yyc[j]/xxc[i])
           zzc_grid[i,j,k]=zzc[k]
        endfor
     endfor
  endfor
  coords={dx:dx,dy:dy,dz:dz,drmin:drmin,dv:dv,ng:ng,nx:nx, ny:ny, nz:nz, $
          xx:xx,yy:yy,zz:zz,xxc:xxc,yyc:yyc,zzc:zzc,$
          rrc_grid:rrc_grid,phi_grid:phi_grid,zzc_grid:zzc_grid, $
          isource:inputs.isource}
  nbi={vinj:vinj,xyz_src:xyz_src,Arot:Arot,Brot:Brot,Crot:Crot}
  ;now make nbigeom for compatibility. write only data for isource,
  ;rest is zero.
  xyz_src2 = fltarr(8,3)
  Arot2 = fltarr(8,3,3)
  Brot2 = fltarr(8,3,3)
  Crot2 = fltarr(8,3,3)
  xyz_src2[isource,*] = xyz_src
  Arot2[isource,*,*] = Arot
  Brot2[isource,*,*] = Brot
  Crot2[isource,*,*] = Crot
  nbgeom={  Arot:   Arot2 $
             , Brot:   Brot2 $
             , Crot:   Crot2 $
             ;; parameters
             , focy: double(focy)     , focz: double(focz) $
             , divy: double(divy)      , divz: double(divz) $
             , bmwidra:double(bmwidra)  , bmwidza:double(bmwidza)  $
             , xyz_src: xyz_src2  }
    


  ;; -------------------------------------------------------
  ;; ------------- LOAD THE SPECTRA ------------------------
  ;; -------------------------------------------------------
  if npa ne 1 and calc_spec eq 1 then begin
     spec_file=result_dir+'/nbi_halo_spectra.bin' 
     openr, 55, spec_file
     readu,55,idum  ;; shot
     readu,55,fdum  ;; time
     readu,55,idum & nlos=idum 
     readu,55,idum & nlambda=idum
     lambda=fltarr(nlambda)
     readu,55,lambda
     fspectra=fltarr(nlambda,nlos)
     readu,55,fspectra
     hspectra=fltarr(nlambda,nlos)
     readu,55,hspectra
     tspectra=fltarr(nlambda,nlos)
     readu,55,tspectra
     halospectra=fltarr(nlambda,nlos)
     readu,55,halospectra
     bremspectra=fltarr(nlambda,nlos)
     readu,55,bremspectra
     close, 55
     ;; ----- NOW LOAD THE FIDA SPECTRUM IF AVAILABLE ------!!
     spec_file=result_dir+'/fida_spectra.bin' 
     if nr_fast gt 10 and file_test(spec_file) then begin
        openr, 55, spec_file
        readu,55,idum  ;; shot
        readu,55,fdum  ;; time
        readu,55,idum & nlos=idum 
        readu,55,idum & nlambda=idum
        lambda_fida=fltarr(nlambda)
        readu,55,lambda
        fidaspectra =fltarr(nlambda,nlos) 
        readu,55,fidaspectra
        close, 55
     endif else fidaspectra=fspectra*0.d0
     spec={ shot:shot, time:time,nlos:nlos $
            ,nlambda:nlambda,lambda: lambda,disp:lambda[1]-lambda[0] $
            ,full: fspectra,half:hspectra $
            ,third:tspectra,halo:halospectra $
            ,fida:fidaspectra,brems:bremspectra }
  endif else spec=0.d0
  
 
  ;; -------------------------------------------------------
  ;; --------- LOAD THE NEUTRAL DENSITY PROFILES -----------
  ;; -------------------------------------------------------
  only_neut:
  neut_file=result_dir+'/neutrals.bin' 
  if file_test(neut_file) then begin
     openr, 55,  neut_file
     readu,55,dum               ;& shot      = fix(dum)     &  print, 'shot  :' ,shot   
     readu,55,dum               ;& time=dum      ;& print, 'time',time
     readu,55,dum  & nx        = fix(dum)  
     readu,55,dum  & ny        = fix(dum)   
     readu,55,dum  & nz        = fix(dum)    
     readu,55,dum  & nlevs     = fix(dum) 
     ;; create arrays
     fdens   =fltarr(nx,ny,nz,nlevs)
     hdens   =fltarr(nx,ny,nz,nlevs)
     tdens   =fltarr(nx,ny,nz,nlevs)
     halodens=fltarr(nx,ny,nz,nlevs)
     ;; read arrays
     readu,55,fdens
     readu,55,hdens
     readu,55,tdens
     readu,55,halodens
     close,55
     neutrals={fdens:fdens, hdens:hdens,tdens:tdens, halodens:halodens,nlevs:nlevs}
     ;; read los information
     if keyword_set(only_neutrals) then begin
        results={neutrals:neutrals}
        return
     endif
  endif else neutrals=0.0
  
  ;; -------------------------------------------------------
  ;; --------------- LOAD THE LOS GEOMETRY -----------------
  ;; -------------------------------------------------------
  los_file=result_dir+'/los.bin' 
  if (calc_spec or calc_wght or npa) and file_test(los_file) then begin
     xyzhead=fltarr(3)
     openr, 55,  los_file
     readu,55 , idum & nchan=idum
     xyzlos        = FLTARR(nchan,3)
     xyzhead       = FLTARR(nchan,3)
     headsize      = FLTARR(nchan) 
     opening_angle = FLTARR(nchan) 
     sigma_pi      = FLTARR(nchan) 
     for i=0,nchan-1 do begin
        readu,55, ddum & xyzhead[i,0]=ddum
        readu,55, ddum & xyzhead[i,1]=ddum
        readu,55, ddum & xyzhead[i,2]=ddum
        readu,55, ddum & xyzlos[i,0]=ddum
        readu,55, ddum & xyzlos[i,1]=ddum
        readu,55, ddum & xyzlos[i,2]=ddum
        readu,55, ddum & headsize[i]=ddum
        readu,55, ddum & opening_angle[i]=ddum
        readu,55, ddum & sigma_pi[i]=ddum
     endfor
     weight=dblarr(nx,ny,nz,nchan)
     readu,55,weight
     close,55
     
     det = set_det( NCHAN=nchan $
                    ,XLOS=reform(xyzlos[*,0]), YLOS=reform(xyzlos[*,1]) $
                    ,ZLOS=reform(xyzlos[*,2]) $
                    ,XHEAD=xyzhead[*,0], YHEAD=xyzhead[*,1] $
                    ,ZHEAD=xyzhead[*,2], HEADSIZE=headsize $
                    ,OPENING_ANGLE=opening_angle, SIGMA_PI=sigma_pi )
     
     detector={det:det,weight:weight}
     rlos=sqrt(xyzlos[*,0]^2+xyzlos[*,1]^2)
     los={nchan:nchan, xyzhead:xyzhead, xyzlos:xyzlos,rlos:rlos $
          ,los_weight:weight, headsize:headsize}
  endif else begin
     det = set_det()
     detector={det:det,weight:replicate(0.d0,coords.nx,coords.ny,coords.nz,1) }
     los={nchan:0}
  endelse



  

  ;; ------------------------------------------------------------
  ;; --------------- LOAD THE PLASMA PARAMETERS -----------------
  ;; ------------------------------------------------------------
  plasma_file=result_dir+'/plasma.bin' 
  if file_test(plasma_file) then begin
     bx=dblarr(nx,ny,nz)
     by=dblarr(nx,ny,nz)
     bz=dblarr(nx,ny,nz)
     b=dblarr(nx,ny,nz,3)
     ex=dblarr(nx,ny,nz)
     ey=dblarr(nx,ny,nz)
     ez=dblarr(nx,ny,nz)
     e=dblarr(nx,ny,nz,3)
     vti=dblarr(nx,ny,nz)
     dene=dblarr(nx,ny,nz)
     denp=dblarr(nx,ny,nz)
     denf=dblarr(nx,ny,nz)
     deni=dblarr(nx,ny,nz)
     te=dblarr(nx,ny,nz)
     ti=dblarr(nx,ny,nz)
     vrotx=dblarr(nx,ny,nz)
     vroty=dblarr(nx,ny,nz)
     vrotz=dblarr(nx,ny,nz)
     vrot=dblarr(nx,ny,nz,3)
     rho_grid=dblarr(nx,ny,nz)
     zeff=dblarr(nx,ny,nz)
     nx=0L
     ny=0L
     nz=0L
     ;; open and read the file
     openr, 55, plasma_file
     readu,55 ,nx
     readu,55 ,ny
     readu,55 ,nz
     readu,55 ,te
     readu,55 ,ti
     readu,55 ,dene
     readu,55 ,denp
     readu,55 ,deni
     readu,55 ,denf
     readu,55 ,vrotx & vrot[*,*,*,0]=vrotx
     readu,55 ,vroty & vrot[*,*,*,1]=vroty
     readu,55 ,vrotz & vrot[*,*,*,2]=vrotz
     readu,55 ,zeff
     readu,55 ,bx & b[*,*,*,0]=bx
     readu,55 ,by & b[*,*,*,1]=by
     readu,55 ,bz & b[*,*,*,2]=bz
     readu,55 ,ex & e[*,*,*,0]=ex
     readu,55 ,ey & e[*,*,*,1]=ey
     readu,55 ,ez & e[*,*,*,2]=ez
     readu,55 ,rho_grid
     close,55
     vtor = sqrt(total(vrot^2,4)) ;vtor = abs. value,  vrot = vector

     plasma={rho_grid:rho_grid,te:te,ti:ti,vtor:vtor,vrot:vrot $
             ,dene:dene,deni:deni,denp:denp,denf:denf $
             ,vti:vti, B:B, E:E,zeff:zeff}
  endif else plasma=0.

  

  ;; -----------------------------------------------------------------
  ;; ----------------- LOAD THE NPA DATA -----------------------------
  ;; -----------------------------------------------------------------
  if npa eq 1 then begin
     print, 'read NPA'
     npa_file =result_dir+'/npa.bin' 
     openr, 55, npa_file
     readu,55,idum             ;& print, 'shot:', fdum
     readu,55,fdum             ;& print, 'time:', fdum
     readu,55,idum & nr_npa=idum
     readu,55,idum & counter=idum

     help,nr_npa,counter

     print, 'nr of particles at detector: ', counter
     if counter gt 0 then begin
        npaipos=fltarr(counter,3)
        npafpos=fltarr(counter,3)
        npav  = fltarr(counter,3)
        npawght=fltarr(counter)
        readu,55,npaipos
        readu,55,npafpos
        readu,55,npav
        readu,55,npawght 
        close,55
        npa={npaipos:npaipos,npafpos:npafpos,npav:npav, npawght:npawght}
     endif else begin
        print, 'attention, no particle reached the detector!'
        stop
     endelse
  endif else npa={npaipos:0}


  ;;---------------------------------------------------
  ;; READ TRANSP FBM FILE STORED IN RESULTS DIRECTORY
  ;;---------------------------------------------------
  if not keyword_set(no_fbm) then begin
     read_fbm,result_dir,fbm_struct,pitch_sign_convention 
  endif else begin
     fbm_struct = {not_loaded:0}
  endelse

  ;; -----------------------------------------------------------------
  ;; ----------------- LOAD THE WEIGHT FUNCTION ----------------------
  ;; -----------------------------------------------------------------
  file=result_dir+'/weight_function.bin' 
  if calc_wght and file_test(file) then begin
  ;; Open file for the outputs
    openr, 55,  file
    readu,55 , idum & shot       = idum
    readu,55 , fdum & time       = fdum 
    readu,55 , idum & ichan_wght = idum
    readu,55 , idum & nen        = idum
    readu,55 , fdum & dE     = fdum
    readu,55 , idum & npitch = idum
    readu,55 , fdum & dpitch = fdum
    readu,55 , idum & nwav   = idum
    readu,55 , fdum & dwav   = fdum
    readu,55 , fdum & wavel_start = fdum
    readu,55 , fdum & wavel_end   = fdum
    readu,55 , idum & nchan       = idum   
    angle=fltarr(nchan)
    radius=fltarr(nchan)
    weight_tot=fltarr(nchan,nwav,nen,npitch)
    dummy_arr=fltarr(nwav,nen,npitch)
    for chan=0,nchan-1 do begin
       readu,55 , idum 
       readu,55 , fdum & angle[chan]=fdum
       readu,55 , fdum & radius[chan]=fdum
       readu,55, dummy_arr
       weight_tot[chan,*,*,*]=dummy_arr
    endfor
    close,55
    central_wavel=(findgen(nwav)+0.5)*dwav+wavel_start
    energyarr=(findgen(nen)+0.5)*dE
    pitcharr=(findgen(npitch)+0.5)*dpitch-1.
    ;in fidasim binary file, pitch is defined in comp. to B-field.
    ;now revert to more common convention (to current) -> - sign.
    pitcharr*=pitch_sign_convention
    emax=max(energyarr)+0.5*dE

  	shot=inputs.shot
  	time=inputs.time
  	if inputs.calc_wght eq 1 then begin
     	wfunct=weight_tot
     	fbm=fbm_struct
     	nen_transp=n_elements(fbm.energy)
     	npitch_transp=n_elements(fbm.pitch) 
     	dE_transp=fbm.energy[1]-fbm.energy[0]
     	dpitch_transp=abs(fbm.pitch[1]-fbm.pitch[0])
     	radiance= fltarr(nchan,nwav)
     	fbm_tot = fltarr(nchan,nwav)
     	mean_fbm= fltarr(nen_transp,npitch_transp,nchan)
  	endif    ;;Calculate synthetic spectra

    for ichan=0, los.nchan-1 do begin
     rlos=sqrt(los.xyzlos[ichan,0]^2 + $
               los.xyzlos[ichan,1]^2)/100. ;m
     zlos=los.xyzlos[ichan,2]/100. ;[m]
     if inputs.calc_wght eq 1 then begin
        if ichan_wght ne ichan and $
           ichan_wght gt 0 then continue
        ;; ------------------------------------------
        ;;------------ CALCUATE mean value of fbm----
        ;; ------------------------------------------
        rad=0.d0
        for i=0,coords.nx-1 do begin
           for j=0,coords.ny-1 do begin
              for k=0,coords.nz-1 do begin
                 los_wght=detector.weight[i,j,k,ichan]
                 if los_wght gt 0. then begin
                    ;; determine mean values like the halo density along LOS
                    wght=(  neutrals.fdens[i,j,k,2] + $
                            neutrals.hdens[i,j,k,2] + $
                            neutrals.tdens[i,j,k,2] + $
                            neutrals.halodens[i,j,k,2]) * los_wght
                    rad=rad+wght
                    
                    rrc=sqrt(coords.xxc[i]^2+coords.yyc[j]^2)
                    zzc=coords.zzc[k]
                    dr=2.       ;[cm]
                    dz=2.       ;[cm]
                    dummy=min((fbm.r2d-rrc+dr)^2+(fbm.z2d-zzc+dz)^2,fbm_index1)
                    dummy=min((fbm.r2d-rrc-dr)^2+(fbm.z2d-zzc+dz)^2,fbm_index2)
                    dummy=min((fbm.r2d-rrc+dr)^2+(fbm.z2d-zzc-dz)^2,fbm_index3)
                    dummy=min((fbm.r2d-rrc-dr)^2+(fbm.z2d-zzc-dz)^2,fbm_index4)
                    dummy=min((fbm.r2d-rrc)^2+(fbm.z2d-zzc)^2,fbm_index5)
                    mean_fbm[*,*,ichan] = mean_fbm[*,*,ichan] + ( $
                                       fbm.fbm[*,*,fbm_index1] + $
                                       fbm.fbm[*,*,fbm_index2] + $
                                       fbm.fbm[*,*,fbm_index3] + $
                                       fbm.fbm[*,*,fbm_index4] + $
                                       fbm.fbm[*,*,fbm_index5])/5.*wght 
                 endif
              endfor
           endfor
        endfor
        mean_fbm[*,*,ichan]=mean_fbm[*,*,ichan]/rad
        ;;------------------------------------------------------------
        ;; map FBM on the energy and pitch grid of the weight function
        ;;------------------------------------------------------------
        mean_fbm2=fltarr(nen,npitch)
        for ie=0,nen-1 do begin
           dummy=min(abs(fbm.energy-energyarr[ie]),eindex)
           for ip=0,npitch-1 do begin
              dummy=min(abs(fbm.pitch-pitcharr[ip]),pindex)
              mean_fbm2[ie,ip]=mean_fbm[eindex,pindex,ichan]
              ;;[fast_ion/cm^3/dP/dE]
           endfor
        endfor
        ;;------------------------------------------------------------
        ;; ------ CALCULATE SYNTHETIC SPECTRA AND PROFIELES ----------
        ;;------------------------------------------------------------
        for ii=0,nwav-1 do begin
           ;; PRODUCT with fast-ion distribution funciton
           prod=replicate(0.,nen,npitch)
           for ie=0,nen-1 do begin
              for ip=0,npitch-1 do begin
                 prod[ie,ip]=mean_fbm2[ie,ip]*wfunct[ichan,ii,ie,ip]
                 ;;--> [ph/(s cm^2 dP keV)]
              endfor
           endfor
           radiance[ichan,ii]=total(prod[*,*])*dE * $
                              dpitch/(4.d0*!pi)*1.d4 / $
                              dwav ;;--> [ph/(s m^2 sr nm)]
        endfor
     endif
  endfor ;; loop over channels

    wfunct={nchan:nchan,ichan_wght:ichan_wght,nen:nen $
            ,dE:dE,emax:emax,emin:0.,npitch:npitch,dpitch:dpitch   $
            ,nwav:nwav,dwav:dwav,wavel_start:wavel_start,wavel_end:wavel_end   $
            ,central_wavel:central_wavel,wght_spec:transpose(radiance),energyarr:energyarr,pitcharr:pitcharr $
            ,weight_tot:weight_tot   $
            ,angle:angle,radius:radius}
 endif else wfunct={ichan_wght:0}


  ;; birth profile
  birth_file=result_dir+'/birth.bin' 
  if file_test(birth_file) then begin
     openr, 55,  birth_file
     readu,55,fdum & shot=fdum
     readu,55,fdum & time=fdum
     readu,55,fdum & nx=fdum
     readu,55,fdum & ny=fdum
     readu,55,fdum & nz=fdum
     readu,55,fdum & npitch=fdum
     ;; create arrays
     birth_dens=fltarr(nx,ny,nz,3,npitch)
     readu,55,birth_dens
     close,55
  endif else birth_dens=0.d0  


;store everything in result structure:
  results={inputs:inputs,coords:coords $
           ,nbgeom:nbgeom,detector:detector,plasma:plasma $
           ,neutrals:neutrals,npa:npa,spec:spec $
           ,wfunct:wfunct, FBM:FBM_struct,birth_dens:birth_dens $
           ,nbi:nbi,los:los  }
end








    

