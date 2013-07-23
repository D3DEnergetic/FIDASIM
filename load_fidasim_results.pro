pro load_fidasim_results,results,result_dir,only_inputs=only_inputs
  ;; routine to load the results from fidasim2.0
  ;; read inputs
  print,'load fidasim results'
  input_file=result_dir+'/inputs.dat' 
  
  sdum=''
  idum=1
  ldum=1L
  fdum=1.d0
  transp_runid=''
  fidasim_runid=''
  diag=''
  openr,55,input_file
  readf,55,sdum                 ;& print, sdum
  readf,55,sdum                 ;& print, sdum
  readf,55,ldum & shot=ldum & print,'shot: ',  shot
  readf,55,fdum & time=float(fdum) & print,'time: ', time
  readf,55,fidasim_runid
  readf,55,diag 
  diag=strmid(diag,1,3)
  readf,55,sdum                 ;& print, sdum
  readf,55,idum & nospec=idum
  readf,55,idum & nofida=idum               
  readf,55,idum & npa=idum      ;& print, 'npa: ', npa
  readf,55,idum & load_neutrals=idum
  readf,55,idum & guidingcenter=idum
  readf,55,idum & f90brems=idum
  readf,55,idum & calc_wght=idum
  readf,55,sdum;'# weight function settings:'
  readf,55,sdum;inputs.nr_wght,f='(i9,"      # number velocities")'
  readf,55,sdum;inputs.ichan_wght,f='(i3,"      # channel for weight function")'
  readf,55,sdum;inputs.emax_wght,f='(1f12.2,"       # emax for weights")'
  readf,55,sdum;inputs.dwav_wght,f='(1f12.5,"       # dwav")'
  readf,55,sdum;inputs.wavel_start_wght,f='(1f12.5,"       # wavel_start")'
  readf,55,sdum;inputs.wavel_end_wght,f='(1f12.5,"       # wavel_end")'
  readf,55,sdum;'# Monte Carlo settings:'
  readf,55,fdum & nr_fida=long(fdum)
  readf,55,fdum & nr_ndmc=long(fdum)
  readf,55,fdum & nr_halo=long(fdum)
  readf,55,sdum;impurity_charge,f='(i2,"             # Impurity charge")'
  readf,55,sdum;'# Location of transp cdf file:'
  readf,55,transp_runid 
  readf,55,sdum                 ;& print, sdum
  readf,55,sdum                 ;& print, sdum
  readf,55,fdum & ai=fdum       ;& print, 'ai: ', ai
  readf,55,fdum & ab=fdum       ;& print, 'ab: ', ai
  for i=0,4 do readf,55,sdum
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
  for i=0,2 do readf,55,sdum
  nr_src=1
  xyz_src=dblarr(3)
  Arot   =dblarr(3,3)
  Brot   =dblarr(3,3)
  Crot   =dblarr(3,3)
  for i=0,8 do readf,55,sdum
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
  
  inputs={shot: shot, time: time,diag:diag $
          , transp_runid:transp_runid $
          , fidasim_runid:fidasim_runid $
          , nr_fida:nr_fida,nr_ndmc:nr_ndmc,nr_halo:nr_halo $
          , nofida:nofida,load_neutrals:load_neutrals,guidingcenter:guidingcenter,f90brems:f90brems,calc_wght:calc_wght $
          , nx:nx, ny:ny, nz:nz $
          , x0:xx[0],x1:xx[nx-1]+dx $
          , y0:yy[0],y1:yy[ny-1]+dy $
          , z0:zz[0],z1:zz[nz-1]+dz}
  if keyword_set(only_inputs) then begin
     results=inputs
     return
  endif
  coords={nx:nx, ny:ny, nz:nz,dx:dz,dy:dy,dz:dz,xx:xx,yy:yy,zz:zz $
          ,xxc:xx+.5*dx,yyc:yy+.5*dy,zzc:zz+.5*dz}
  nbi={vinj:vinj,xyz_src:xyz_src,Arot:Arot,Brot:Brot,Crot:Crot}
  


  if npa ne 1 and nospec eq 0 then begin
     ;; load spectra
     spec_file=result_dir+'/nbi_halo_spectra.dat' 
     openr, 55, spec_file
     readf,55,idum  ;; shot
     readf,55,fdum  ;; time
     readf,55,diag  ;; diag
     readf,55,idum  & nlos      = fix(idum) 
     readf,55,idum  & nlambda   = fix(idum) 
     readf,55,sdum
     lambda=fltarr(nlambda)
     for i=0,nlambda-1 do begin
        readf,55,fdum & lambda[i]=fdum
     endfor
     readf,55,sdum 
     fspectra=fltarr(nlambda,nlos)
     hspectra=fltarr(nlambda,nlos)
     tspectra=fltarr(nlambda,nlos) 
     halospectra =fltarr(nlambda,nlos) 
     bremspectra =fltarr(nlambda,nlos) 
     for i=0,nlambda-1 do begin
        for los=0,nlos-1 do begin
           readf,55,fdum & fspectra[i,los]   = fdum
           readf,55,fdum & hspectra[i,los]   = fdum
           readf,55,fdum & tspectra[i,los]   = fdum
           readf,55,fdum & halospectra[i,los]= fdum
           readf,55,fdum & bremspectra[i,los]= fdum
        endfor
     endfor        
     close, 55
     
     if nofida ne 1 then begin
        spec_file=result_dir+'/fida_spectra.dat' 
        openr, 55, spec_file
        readf,55,idum 
        readf,55,fdum  
        readf,55,sdum
        readf,55,idum  
        readf,55,idum  
        readf,55,sdum
        lambda_fida=fltarr(nlambda)
        for i=0,nlambda-1 do begin
           readf,55,fdum & lambda[i]=fdum
        endfor
        readf,55,sdum 
        fidaspectra =fltarr(nlambda,nlos) 
        for i=0,nlambda-1 do begin
           for los=0,nlos-1 do begin
              readf,55,fdum & fidaspectra[i,los]  = fdum  
           endfor
        endfor   
        close, 55
     endif else fidaspectra=fspectra*0.d0
     
     spec={ shot:shot, time:time,lambda: lambda,disp:lambda[1]-lambda[0] $
            ,full: fspectra,half:hspectra $
            ,third:tspectra,halo:halospectra $
            ,fida:fidaspectra,brems:bremspectra }
  endif else spec=0.d0
  
  
  ;; read neutral density
  neut_file=result_dir+'/neutrals.bin' 
  openr, 55,  neut_file
  readu,55,dum         ;& shot      = fix(dum)     &  print, 'shot  :' ,shot   
  readu,55,dum         ;& time=dum      ;& print, 'time',time
  readu,55,dum         ;& nx        = fix(dum)  
  readu,55,dum         ;& ny        = fix(dum)   
  readu,55,dum         ;& nz        = fix(dum)    
  readu,55,dum  & nlevs     = fix(dum) 

  ;; create arrays
  fdens=fltarr(nx,ny,nz,nlevs)
  hdens=fltarr(nx,ny,nz,nlevs)
  tdens=fltarr(nx,ny,nz,nlevs)
  halodens=fltarr(nx,ny,nz,nlevs)
  for ix=0,nx-1 do begin
     for iy=0,ny-1 do begin
        for iz=0,nz-1 do begin
           for jj=0,nlevs-1 do begin
              readu,55, dum & fdens[ix,iy,iz,jj]=dum 
              readu,55, dum & hdens[ix,iy,iz,jj]=dum
              readu,55, dum & tdens[ix,iy,iz,jj]=dum
              readu,55, dum & halodens[ix,iy,iz,jj]=dum   
           endfor                   
        endfor
     endfor
  endfor
  close,55
  neutrals={fdens:fdens, hdens:hdens,tdens:tdens, halodens:halodens,nlevs:nlevs}
  ;; read los information


  if  nospec eq 0 then begin
     nchan=0L
     fdummy=1.e0
     dummy=1.d0
     xyzhead=dblarr(3)
     los_file=result_dir+'/los.bin' 
     openr, 55,  los_file
     readu,55 , nchan
     xyzlos=dblarr(nchan,3)
     xyzhead=dblarr(nchan,3)
     headsize=dblarr(nchan) ;; size of NPA head
     for i=0,nchan-1 do begin
        readu,55 , dummy & xyzhead[i,0]=dummy
        readu,55 , dummy & xyzhead[i,1]=dummy
        readu,55 , dummy & xyzhead[i,2]=dummy
        readu,55 , dummy & headsize[i]=dummy
        readu,55, dummy & xyzlos[i,0]=dummy
        readu,55, dummy & xyzlos[i,1]=dummy
        readu,55, dummy & xyzlos[i,2]=dummy
     endfor
     readu,55 , dummy
     los_weight=fltarr(nx,ny,nz,nchan)
     for ix=0,nx-1 do begin
        for iy=0,ny-1 do begin
           for iz=0,nz-1 do begin       
              for ilos=0,nchan-1 do begin
                 readu ,55 , fdummy
                 los_weight[ix,iy,iz,ilos]=fdummy
              endfor  
           endfor
        endfor
     endfor
     close,55
     los={nchan:nchan, xyzhead:xyzhead, xyzlos:xyzlos,los_weight:los_weight}
  endif else los={nchan:0, xyzhead:0, xyzlos:0,los_weight:0}
  

  
  ;;read plasma
  B=fltarr(nx,ny,nz,3)
  E=fltarr(nx,ny,nz,3)
  vti=fltarr(nx,ny,nz)
  dene=fltarr(nx,ny,nz)
  denp=fltarr(nx,ny,nz)
  denf=fltarr(nx,ny,nz)
  deni=fltarr(nx,ny,nz)
  te=fltarr(nx,ny,nz)
  ti=fltarr(nx,ny,nz)
  vrot=fltarr(nx,ny,nz,3)
  rho=fltarr(nx,ny,nz)
  nx=0L
  ny=0L
  nz=0L

  ;; define dummy variables
  dummy=double(1.d0)
  dumB1 =double(1.d0)
  dumB2 =double(1.d0)
  dumB3 =double(1.d0)
  dumE1 =double(1.d0)
  dumE2 =double(1.d0)
  dumE3 =double(1.d0)
  dumvti =double(1.d0)
  dumte =double(1.d0)
  dumti =double(1.d0)
  dumdene =double(1.d0)
  dumdenp =double(1.d0)
  dumdeni =double(1.d0)
  dumdenf =double(1.d0)
  dumzeff =double(1.d0)
  dumvrot1 =double(1.d0)
  dumvrot2 =double(1.d0)
  dumvrot3 =double(1.d0)

  dumrho =double(1.d0)
  ;; open and read the file
  plasma_file=result_dir+'/plasma.bin' 
  openr, 55, plasma_file
  readu,55 , nx
  readu,55 , ny
  readu,55 , nz
  for ix=0,nx-1 do begin
     for iy=0,ny-1 do begin
        for iz=0,nz-1 do begin
           readu,55 , dumte  , dumti    $
                 , dumdene  , dumdenp  $
                 , dumdeni  , dumvrot1  $
                 , dumvrot2  , dumvrot3  $
                 , dumB1   , dumB2   $
                 , dumB3   , dumE1  $
                 , dumE2  , dumE3  $
                 , dumrho   , dumdenf $
                 , dumzeff
           vti[ix,iy,iz]= sqrt(dumti*1.d3*e0/(ai*mass_u))*1.d2 ; [cm/s]
           B[ix,iy,iz,*]=[dumB1,dumB2,dumB3]
           E[ix,iy,iz,*]=[dumE1,dumE2,dumE3]
           te[ix,iy,iz]=dumte
           ti[ix,iy,iz]=dumti    
           dene[ix,iy,iz]=dumdene
           denf[ix,iy,iz]=dumdenf
           denp[ix,iy,iz]=dumdenp
           deni[ix,iy,iz]=dumdeni
           vrot[ix,iy,iz,*]=[dumvrot1,dumvrot2,dumvrot3]
           rho[ix,iy,iz]=dumrho
        endfor
     endfor
  endfor
  close,55
  plasma={te:te,ti:ti,vrot:vrot,dene:dene,deni:deni,denp:denp,denf:denf $
          ,vti:vti, rho:rho,B:B, E:E}
  
  if npa eq 1 then begin
     print, 'read NPA'
     npa_file =result_dir+'/npa.bin' 
     dummy=float(1.0)
     openr, 55, npa_file
     readu,55,dummy             ;& print, 'shot:', dummy
     readu,55,dummy             ;& print, 'time:', dummy
     readu,55,dummy & counter=long(dummy)
     npaipos=fltarr(counter,3)
     npafpos=fltarr(counter,3)
     npav  = fltarr(counter,3)
     npawght=fltarr(counter)
     for i =0L, counter-1 do begin
        readu,55,dummy & npaipos[i,0]=dummy
        readu,55,dummy & npaipos[i,1]=dummy
        readu,55,dummy & npaipos[i,2]=dummy
        readu,55,dummy & npafpos[i,0]=dummy
        readu,55,dummy & npafpos[i,1]=dummy
        readu,55,dummy & npafpos[i,2]=dummy
        readu,55,dummy & npav[i,0]  =dummy
        readu,55,dummy & npav[i,1]  =dummy
        readu,55,dummy & npav[i,2]  =dummy
       readu,55,dummy & npawght[i] =dummy
    endfor
     close,55
     npa={npaipos:npaipos,npafpos:npafpos,npav:npav, npawght:npawght}
  endif else npa={npaipos:0}

  if calc_wght eq 1 then begin
  ;; Open file for the outputs
    file=result_dir+'/weight_function.bin' 
    fdum=1.e0
    dummy=1.d0
    openr, 55,  file
    readu,55 , fdum & shot=fdum ;& print ,shot
    readu,55 , fdum & time=fdum ;& print, time
    readu,55 , fdum & ichan=fdum-1 & print,'chan selection:',ichan
    readu,55 , fdum & nen=fdum  ;& print, nen
    readu,55 , fdum & dE=fdum
    readu,55 , fdum & npitch=fdum
    readu,55 , fdum & dpitch=fdum
    readu,55 , fdum & nwav=fdum
    readu,55 , fdum & dwav=fdum
    readu,55 , fdum & wavel_start=fdum
    readu,55 , fdum & wavel_end=fdum
    readu,55 , fdum & nchan=fdum   
    angle=fltarr(nchan)
    radius=fltarr(nchan)
    weight_tot=fltarr(nchan,nwav,nen,npitch)
    for chan=0,nchan-1 do begin
       readu,55 , fdum & print, 'ichan: ',long(fdum)
       readu,55 , fdum & angle[chan]=fdum
       readu,55 , fdum & radius[chan]=fdum
       for ii=0,nwav-1 do begin   
          for i=0,nen-1 do begin
             for j=0,npitch-1 do begin   
                readu,55 , fdum & weight_tot[chan,ii,i,j]=fdum
             endfor
          endfor
       endfor
    endfor
    close,55
    central_wavel=(findgen(nwav)+0.5)*dwav+wavel_start
    energyarr=(findgen(nen)+0.5)*dE
    pitcharr=(findgen(npitch)+0.5)*dpitch-1.
    emax=max(energyarr)
    wfunct={nchan:nchan,ichan:ichan,nen:nen $
            ,dE:dE,emax:emax,npitch:npitch,dpitch:dpitch   $
            ,nwav:nwav,dwav:dwav,wavel_start:wavel_start,wavel_end:wavel_end   $
            ,central_wavel:central_wavel,energyarr:energyarr,pitcharr:pitcharr $
            ,weight_tot:weight_tot   $
            ,angle:angle,radius:radius}

 endif else wfunct={ichan:0}

  
  results={inputs:inputs,coords:coords $
           ,nbi:nbi,los:los,plasma:plasma $
           ,neutrals:neutrals,npa:npa,spec:spec $
           ,wfunct:wfunct}

  ;save,filename='fidasim_results.idl',results
end





    

