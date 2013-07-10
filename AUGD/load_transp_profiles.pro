pro read_transp,shot,time_in,kind,diag,rho,output,rhostr,err
  directory='/u/bgeiger/tr_client/AUGD/'+string(shot,f='(i5)')+'/'
  file=directory+kind+string(shot,f='(i5)')+'.'+diag

  openr, 55, file ,error = err  ; open ufile
  IF (err ne 0) then return
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
  while cc lt nmax do begin
     READF, 55, line
     dummy = strsplit(line,' ',/extract) 
     ll++
     cc1=cc
     if n_elements(dummy) eq 6 then begin
        data_raw[indgen(6)+cc]=double(dummy)
        cc=cc+6
     endif else begin
        for i=0,n_elements(dummy)-1 do begin
           if strlen(dummy[i]) lt 15 then begin
              data_raw[cc]=dummy[i]
              cc++
           endif else begin
              dummy2 = strsplit(dummy[i],'-',/extract) 
              for ii=0,n_elements(dummy2)-1 do begin 
                 factor=-1.
                 if ii eq 0 then begin
                    mpos=strpos(dummy[i],'-')
                    if mpos ne 0 then factor=1.
                 endif
                 ;; there is a problem when the values become E-01!!
                 if strlen(dummy2[ii]) lt 11 then begin
                    dummy2[ii]=dummy2[ii]+'-'+dummy2[ii+1]
                    ii=ii+1
                 endif
                 data_raw[cc]=factor*dummy2[0] 
                 cc++
              endfor
           endelse
        endfor
     endelse
     
     cc2=cc
     if cc2-cc1 gt 6 then stop
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
  index=where(time_arr gt time_in-0.001 and time_arr lt time_in+0.001,nind)
  if nind gt 1 then begin
     ;; take mean of the data
     output=total(data[index,*],1)/double(nind)
   ;  print, 'time file:', mean(time_arr[index])
  endif else begin
        ;; take the one time point
        dummy=min(abs(time_arr-time_in),index)
        output=reform(data[index,*])
    ;    print, 'time, time file:', time_in, mean(time_arr[index])
    ;    print, 'dt:', time_arr[index]-time_in
  endelse

end




pro load_transp_profiles,shot,time,profiles,rhostr_ne
  doplot=1


 
  if !VERSION.MEMORY_BITS eq 32 then begin 
     defsysv,'!libddww','/afs/ipp/aug/ads/lib/@sys/libddww.so' 
     defsysv,'!libkk'  ,'/afs/ipp/aug/ads/lib/@sys/libkk.so'
  endif else begin
     defsysv,'!libddww','/afs/ipp/aug/ads/lib64/@sys/libddww.so' 
     defsysv,'!libkk'  ,'/afs/ipp/aug/ads/lib64/@sys/libkk.so'
  endelse   


  print, shot, time
  nr_rho=120
  rho=dindgen(nr_rho)/(nr_rho-1.)*1.2
  profiles=    { time:time              $
                 , rho:  rho      $
                 , ti:  dblarr(nr_rho)    $  
                 , vtor: dblarr(nr_rho)   $
                 , te: dblarr(nr_rho)     $
                 , dene: dblarr(nr_rho)   $
                 , zeff:dblarr(nr_rho)   }
  ;; add zero values to rho > 1.!!
  rho_sol=(dindgen(20)+0.5)/20.*0.2+1.
  val_sol=replicate(0.d0,20)
  ;; NE
  err=1
 ; if err ne 0 then read_transp,shot,time,'N','MOD',rho_dene,dene,rhostr_ne,err
  if err ne 0 then read_transp,shot,time,'N','IDA',rho_dene,dene,rhostr_ne,err
  if err ne 0 then read_transp,shot,time,'N','IDZ',rho_dene,dene,rhostr_ne,err
  if err ne 0 then read_transp,shot,time,'N','VTA',rho_dene,dene,rhostr_ne,err
  index=where(dene gt 0.0)   
  rho_dene=[rho_dene[index],rho_sol]
  dene=[dene[index]*1.d6,val_sol]
  profiles.dene =  interpol(dene,rho_dene,profiles.rho)  

  ;; TE
  err=1
  if err ne 0 then  read_transp,shot,time,'E','IDA',rho_te,te,rhostr_te,err
  if err ne 0 then  read_transp,shot,time,'E','IDZ',rho_te,te,rhostr_te,err
  if err ne 0 then  read_transp,shot,time,'E','CEC',rho_te,te,rhostr_te,err
  if err ne 0 then  read_transp,shot,time,'E','VTA',rho_te,te,rhostr_te,err
  index=where(te gt 0.0)  
  rho_te=[rho_te[index],rho_sol]
  te=[te[index],val_sol]     
  profiles.te   =  interpol(te  ,rho_te,profiles.rho)   
 
  if rhostr_ne ne rhostr_te then stop

  ;; ZEFF
  err=1
  if err ne 0 then read_transp,shot,time,'Z','IDZ',rho_zeff,zeff,rhostr_zef,err
  if err ne  0 then begin
     ;; global estimate 
     read_transp,shot,time,'Z','ZEF',rho_zeff,zeff,rhostr_zef,err
     rho_zeff=rho_dene
     if err ne 0 then begin
        zeff=1.5d0
        zeff=replicate(zeff,n_elements(rho_dene))
        rhostr_zef='Zeff   '
     endif
  endif
  index=where(zeff ge 0.0)   
  rho_zeff=[rho_zeff[index],rho_sol]
  zeff=[zeff[index],val_sol]
  profiles.zeff =  interpol(zeff,rho_zeff,profiles.rho)

  if rhostr_zef ne rhostr_te and rhostr_zef ne 'Zeff   ' then stop

  ;; TI
  err=1
  if err ne 0 then read_transp,shot,time,'I','CEZ',rho_ti,ti,rhostr_ti,err
  if err ne 0 then read_transp,shot,time,'I','CHZ',rho_ti,ti,rhostr_ti,err
  index=where(ti gt 0.0)     
  rho_ti=[rho_ti[index],rho_sol]
  ti=[ti[index],val_sol]
  profiles.ti   =  interpol(ti  ,rho_ti,profiles.rho) 
 
  if rhostr_ti ne rhostr_te then stop
  ;;VTOR
  err=1
  if err ne 0 then  read_transp,shot,time,'V','CEZ',rho_vtor,vtor,rhostr_vtor,err
  if err ne 0 then  read_transp,shot,time,'V','CHZ',rho_vtor,vtor,rhostr_vtor,err
  index=where(vtor gt 0.0)   
  rho_vtor=[rho_vtor[index],rho_sol]
  vtor=[vtor[index]*1.d-2,val_sol]
  profiles.vtor   =  interpol(vtor  ,rho_vtor,profiles.rho) 

 if rhostr_vtor ne rhostr_te then stop


  ;; correction if interpolation goes to infinity
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


end








