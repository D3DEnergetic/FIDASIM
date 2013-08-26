@~/FIDASIM/AUGD/rz2rho.pro
@~/FIDASIM/LIB/nstring.pro
pro read_transp,shot,time_in,file,rho,output,rhostr,err $
                ,quantity=quantity, unit=unit, average=average, outtime=outtime
  openr, 55, file ,error = err  ; open ufile
  IF (err ne 0) then begin
     print, 'Warning: transp file not found: ', file
     print, 'Trying other diagnostics...'
     return
     ;;stop
  endif
  line    = ''
  READF, 55, line ;& print, line
  dummy=strsplit(line,' ',/extract)
  shot = (strsplit(dummy[0],'AUGD',/extract))[0]  
  dims=float(dummy[1])
  READF, 55, line ;& print, line 
  READF, 55, line ;& print, line
  READF, 55, line ;& print, line
  READF, 55, line ;& print, line
  rhostr=strmid(line,1,7)
  READF, 55, line ;& print, line
  quantity_unit =line
  quantity=strtrim(strmid(line,1,7),2)
  unit=strtrim(strmid(line,21,10),2)
  if dims eq 2 then begin
     READF, 55, line ;& print, line 
  endif
  ntimes= 1L
  nrho  =1L
  READF, 55, ntimes ;& print, 'ntimes: ', ntimes
  if dims eq 2 then READF, 55, nrho  ;& print, 'nr rho:', nrho
  ;; READ the TIME ARRAY
  time_arr=dblarr(ntimes)
  readf,55,time_arr
  ;; READ the rho ARRAY
  if dims eq 2 then begin
     rho=dblarr(nrho)
     readf,55,rho
  endif 
 
  ll=0L
  nmax=long(ntimes*nrho)
  data_raw=dblarr(nmax)
  cc=0L

  ;new parse-code (hopefully strlen and stringstart are always the same)
  strlen=13  ;length of one double string (including sign)
  stringstart = indgen(6)*strlen+1 

  while cc lt nmax do begin
     READF, 55, line
     dummy = strmid(line, stringstart, strlen)  ;cut line string
     dummy = dummy[where(dummy ne '')]   ;remove empty strings

     if n_elements(dummy) eq 6 then begin
        data_raw[indgen(6)+cc]=double(dummy)
        cc=cc+6
     endif else begin
        for i=0,n_elements(dummy)-1 do begin
           data_raw[cc]=dummy[i]
           cc++
        endfor
     endelse
  endwhile
  close,55 

  
  data=fltarr(ntimes,nrho)
  cc=0L
  for i=0,nrho-1 do begin
     for j=0,ntimes-1 do begin     
        data[j,i]=data_raw[cc]
        cc++
     endfor
  endfor

  
  if n_elements(average) ne 0 then begin
     delta_t = average/2.
  endif else begin
     delta_t=0.001
  endelse


  index=where(time_arr gt time_in-delta_t and time_arr lt time_in+delta_t,nind)
  if nind gt 1 then begin
     ;; take mean of the data
     print, 'averaged over ' + nstring(nind) + ' time points.'
     output=total(data[index,*],1)/double(nind)
   ;  print, 'time file:', mean(time_arr[index])
  endif else begin
        ;; take the one time point
        dummy=min(abs(time_arr-time_in),index)
        output=reform(data[index,*])
    ;    print, 'time, time file:', time_in, mean(time_arr[index])
    ;    print, 'dt:', time_arr[index]-time_in
     endelse
  outtime= mean(time_arr[index])
  if err eq 0 then begin
     print, 'Ufile loaded: ',strmid(file,strlen(file)-21,21)
  endif
end







pro load_transp_profiles,shot,time,profiles,cdf_file
  ;IN: shot, time, cdf_file
  ;OUT: profiles


  ti_diag='' & rot_diag='' & te_diag='' & dene_diag='' & zeff_diag=''
  constant_zeff=0

  ;;cdf_file can be the shot (as string with length=5), or the cdffile.
  ;;->figure out u-file directory from that.
  ;;if cdffile is from SSFPQL/TORIC, then find directory differently:

  is_toric = stregex(cdf_file, 'toric', /fold) ne -1
  if is_toric then begin
     icut= (strsplit(cdf_file, '/'))[-1]
     dir = strmid(cdf_file, 0, icut)
     ;dir with leading /    
     transp_runid=''
     shot_str=string(shot,f='(i5)')
  endif else begin
  ;;->figure out u-file directory

  ;; ----------------------------------------------------------------
  ;; - USE the cdf_file name to define where the u-files are stored -
  ;; ----------------------------------------------------------------
  if strlen(cdf_file) gt 5 then begin
     dummy=strsplit(cdf_file,'tr_client',/regex,/extract)
     home_dir=dummy[0]
     dummy=strsplit(dummy[1],'/',/extract)
     shot_str=dummy[1]
     dir=home_dir+'tr_client/'+dummy[0]+'/'+shot_str+'/'
     transp_runid=dummy[2]+'/'+dummy[3]
     if long(shot_str) ne shot then begin
        print, 'ATTENTION, THE KINETIC PROFILES WILL BE LOADED FOR #' +shot_str
        print, 'AND NOT FOR THE SELECTED DISCHARGE (#'+string(shot,f='(i5)')+')'
        print, 'CHANGE THE .CDF FILE IN start_fidasim.pro OR CONTINUE WITH .c!'
        stop
     endif
     ;; ---------------------------------------------------------
     ;; ------ CHECK TRANSP NAMELIST WHICH U-FILES ARE USED -----
     ;; ---------------------------------------------------------
     file= dir+strmid(transp_runid,0,strlen(transp_runid)-9)+'TR.DAT'
     if file_test(file) then begin
        openr, 55, file         ;open namelist
        line=''
        WHILE ~ EOF(55) DO BEGIN
           READF, 55, line
           if (strsplit(line, '!'))[0] eq 1 then continue   ;comment only line
           commands = (strsplit(line, '!', /extra))[0]      ;comments removed
           commands = strcompress(commands, /remove_all) ;remove blank spaces
           ;; TI DIAG
           if stregex(commands, 'extti2', /fold_case) ne -1 then begin
              dum=strsplit(commands,"'",/extract)
              ti_diag=dum[1]
              continue
           endif
           ;; ROTATION DIAG
           if stregex(commands, 'extvp2', /fold_case) ne -1 then begin
              dum=strsplit(commands,"'",/extract)
              rot_diag=dum[1]
              continue
           endif
           ;; ELECTRON DENSITY DIAG
           if stregex(commands, 'extner', /fold_case) ne -1 then begin
              dum=strsplit(commands,"'",/extract)
              dene_diag=dum[1]
              continue
           endif
           ;; ELECTRON TEMPERATURE DIAG
           if stregex(commands, 'extter', /fold_case) ne -1 then begin
              dum=strsplit(commands,"'",/extract)
              te_diag=dum[1]
              continue
           endif
           ;; check namelist if Zeff is set manually
           if stregex(commands, 'xzeffi', /fold_case) ne -1 then begin
                                ;-> Zeile gefunden!
              zeffi =strsplit(commands,'xzeffi',/extra,count=count $
                              ,/reg,/fold_case)
              zeffi =strsplit(zeffi, '=', /extra)
              zeff_val = float(zeffi)
              print, 'Zeff found in namelist! Set to: ' + string(zeff_val)
              constant_zeff =1
           endif
           ;; 2D  ZEFF FILE
           if constant_zeff eq 0 then begin
              if stregex(commands, 'extzf2', /fold_case) ne -1 then begin
                 dum=strsplit(commands,"'",/extract)
                 zeff_diag=dum[1]
                 continue
              endif
           endif
        ENDWHILE
        close, 55
     endif
  endif else begin
     ;; ----------------------------------------------------------------
     ;; -------- DEFINE THE POSITION OF THE U-FILES from the shot ------
     ;; ----------------------------------------------------------------
     dir='~/tr_client/AUGD/'+string(shot,f='(i5)')+'/'
     transp_runid=''
     shot_str=string(shot,f='(i5)')
  endelse
  endelse
  


  ;; ----------------------------------------------------------------
  ;; --------------- DEFINE THE OUTPUT STRUCTURE --------------------
  ;; ----------------------------------------------------------------
  nr_rho=120
  rho=dindgen(nr_rho)/(nr_rho-1.)*1.2
  profiles=    { time:time              $
                 , rho:  rho      $
                 , rho_str:  ''      $
                 , ti:  dblarr(nr_rho)    $  
                 , vtor: dblarr(nr_rho)   $
                 , te: dblarr(nr_rho)     $
                 , dene: dblarr(nr_rho)   $
                 , zeff:dblarr(nr_rho)   }
  ;; add zero values to rho > 1.!!
  rho_sol=(dindgen(20)+0.5)/20.*0.2+1.
  val_sol=replicate(0.d0,20)


  ;; --------------------------------------------------
  ;; ------------- ELECTRON DENSITY -------------------
  ;; --------------------------------------------------
  err=1
  file=dir+'N'+shot_str+'.'
  if err ne 0 then read_transp,shot,time,file+dene_diag,rho_dene,dene,rhostr,err
  if err ne 0 then read_transp,shot,time,file+'IDA',rho_dene,dene,rhostr,err
  if err ne 0 then read_transp,shot,time,file+'IDZ',rho_dene,dene,rhostr,err
  if err ne 0 then read_transp,shot,time,file+'VTA',rho_dene,dene,rhostr,err
  if err ne 0 then read_transp,shot,time,file+'DPR',rho_dene,dene,rhostr,err
  if err ne 0 then begin
     print, 'Error: no transp Data for ne found'
     stop
  endif
  index=where(dene gt 0.0)  
  if index[0] eq -1 then stop
  rho_dene=[rho_dene[index],rho_sol]
  dene=[dene[index]*1.d6,val_sol]
  profiles.rho_str=rhostr
  profiles.dene =  interpol(dene,rho_dene,profiles.rho)  

  ;; --------------------------------------------------
  ;; ----------- ELECTRON TEMPERATUERE ----------------
  ;; --------------------------------------------------
  err=1
  file=dir+'E'+shot_str+'.'
  if err ne 0 then  read_transp,shot,time,file+te_diag,rho_te,te,rhostr,err
  if err ne 0 then  read_transp,shot,time,file+'IDA',rho_te,te,rhostr,err
  if err ne 0 then  read_transp,shot,time,file+'IDZ',rho_te,te,rhostr,err
  if err ne 0 then  read_transp,shot,time,file+'CEC',rho_te,te,rhostr,err
  if err ne 0 then  read_transp,shot,time,file+'VTA',rho_te,te,rhostr,err
  if err ne 0 then begin
     print, 'Error: no transp Data for Te found'
     stop
  endif
  index=where(te gt 0.0)  
  if index[0] eq -1 then stop
  rho_te=[rho_te[index],rho_sol]
  te=[te[index],val_sol]     
  profiles.te   =  interpol(te  ,rho_te,profiles.rho)   
  if rhostr ne profiles.rho_str then stop


  ;; --------------------------------------------------
  ;; --------------- ZEFF -----------------------------
  ;; --------------------------------------------------
  if constant_zeff eq 1 then begin
     ;; constang ZEFF value found in NAMELIST
     rho_zeff=rho_dene
     zeff=replicate(zeff_val, n_elements(rho_zeff)) 
  endif else begin
     err=1
     file=dir+'Z'+shot_str+'.'
     if err ne 0 then read_transp,shot,time,file+zeff_diag,rho_zeff,zeff,rhostr,err
     if err ne 0 then read_transp,shot,time,file+'IDZ',rho_zeff,zeff,rhostr,err
     if err ne  0 then begin  ;; load the 1D ZEFF estimate
       ; read_transp,shot,time,file+'ZEF',rho_zeff,zeff1d,rhostr,err
        if err ne 0 then begin
           print, 'No transp Data for Zeff found, set to 1.5!'
           zeff1d=1.5d0
        endif
        rho_zeff=rho_dene
        zeff=replicate(zeff1d, n_elements(rho_zeff)) 
     endif
  endelse
  ;; now add some scrape-off layer data ponints that are set to zero
  index=where(zeff ge 0.0,nzeff)   
  if nzeff LE 0 then stop
  rho_zeff=[rho_zeff[index],rho_sol]
  zeff=[zeff[index],val_sol]
  profiles.zeff =  interpol(zeff,rho_zeff,profiles.rho) 
  if rhostr ne profiles.rho_str and rhostr ne 'Zeff   ' then stop

  ;; --------------------------------------------------
  ;; --------- ION TEMPERATURE ------------------------
  ;; --------------------------------------------------
  err=1  
  file=dir+'I'+shot_str+'.'
  if err ne 0 then read_transp,shot,time,file+ti_diag,rho_ti,ti,rhostr,err
  if err ne 0 then read_transp,shot,time,file+'CEZ',rho_ti,ti,rhostr,err
  if err ne 0 then read_transp,shot,time,file+'CHZ',rho_ti,ti,rhostr,err
  if err ne 0 then begin
     rho_ti=rho_te
     ti=te    
     print,'NO TI FOUND! Using TI eq TE!'
  endif else begin
     index=where(ti gt 0.0)   
     rho_ti=[rho_ti[index],rho_sol]
     ti=[ti[index],val_sol]
  endelse
  profiles.ti   =  interpol(ti  ,rho_ti,profiles.rho) 
  if rhostr ne profiles.rho_str then stop


  ;; --------------------------------------------------
  ;; --------- PLASMA ROTATION ------------------------
  ;; --------------------------------------------------
  err=1
  file=dir+'V'+shot_str+'.'
  if err ne 0 then read_transp,shot,time,file+rot_diag,rho_vtor,vtor,rhostr,err, unit=unit
  if err ne 0 then read_transp,shot,time,file+'CEZ',rho_vtor,vtor,rhostr,err, unit=unit
 if err ne 0 then read_transp,shot,time,file+'CHZ',rho_vtor,vtor,rhostr,err, unit=unit
  if err ne 0 then begin
     print, 'Error: no transp Data for vtor found'
     rho_vtor=rho_te
     vtor=te*0.d0
     print,'NO VTOR FOUND! Using VTOR=0!'
  endif
  if unit ne '(cm/s)' then begin
     print, 'Error: vtor has unit: ' + unit
     print, 'Expected: (cm/s) !'
     stop
  endif
  index=where(vtor) ;; take all nonzero values (pos. and negative)
  rho_vtor=[rho_vtor[index],rho_sol]
  vtor=[vtor[index]*1.d-2,val_sol]
  profiles.vtor   =  interpol(vtor  ,rho_vtor,profiles.rho) 
  if rhostr ne profiles.rho_str then stop




  ;; -----------------------------------------------------
  ;; ---- correction if interpolation goes to infinity ---
  ;; -----------------------------------------------------
  index=where(finite(profiles.vtor) eq 0,nind)
  if nind gt 0 then begin
      if mean(profiles.rho[index]) lt 0.1 then begin
        profiles.vtor[index]=dindgen(nind)*0.d0+vtor[0]
        profiles.ti[index]=dindgen(nind)*0.d0+ti[0]
        profiles.zeff[index]=dindgen(nind)*0.d0+zeff[0]
        profiles.te[index]=dindgen(nind)*0.d0+te[0]
        profiles.dene[index]=dindgen(nind)*0.d0+dene[0]
     endif else stop
  endif
  save,filename='profiles_28061_1.6s.idl',profiles
end








