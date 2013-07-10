pro nbi_data,inputs,einj,pinj,full,half,third,doplot=doplot,ps=ps
  if keyword_set(doplot) then begin
     if not(keyword_set(ps)) then window,1,colors=15,xsize=520,ysize=900
     !P.multi=[0,0,8]
  endif
   
  time_start=float(inputs.time-0.005)
  time_end  =float(inputs.time+0.005)
  shot=long(inputs.shot)
  exp='AUGD'
  diag='NIS'
  
;=============================
; Get power
;=============================
  tbeg=0.
  tend=10.
  name='PNIQ    '
  n_box1=4L
  n_box2=4L
  n_beams=n_box1+n_box2
  ntim=100000L
  tim=fltarr(ntim)
  data=fltarr(ntim,n_beams)
  ctr=3L
  anzsig=1L
  physdim='            '
  timbas=diag+'.        '
  error=0L
  ed=0L
  len=long(ntim)
  ind=[100L,1L,1L]
  status=call_external(!libddww,'ddgetaug',   $
                       'ddgetsgr',error,exp,diag,shot,ed,name,tbeg,tend, $
                       timbas,tim,data,len,ind,anzsig,physdim)
  status=call_external(!libddww,'ddgetaug',   $
                       'xxerror',error,ctr,' ')
  
  dummy=min(abs(tim-time_start),time_start_frame)
  dummy=min(abs(tim-time_end),time_end_frame)
  dat=data[time_start_frame:time_end_frame,*]
  time_dat=tim[time_start_frame:time_end_frame]
  
  plot_dat=data[where(tim gt time_start-0.5 and tim lt time_end+0.5),*]
  time_plot_dat=tim[where(tim gt time_start-0.5 and tim lt time_end+0.5)]
  
  index=[0,2,4,6,1,3,5,7]
  pinj=dblarr(8)
  for i=0,7 do begin
     j=index[i]
     if keyword_set(doplot) then begin
        plot,time_plot_dat, plot_dat[*,j]
        oplot,time_dat,dat[*,j],color=254
     endif
     
     dummy=where(dat[*,j] gt 5.e5,nelems)
     if nelems gt 0 then begin
        pinj[i]=mean(dat(dummy,j))/1.d6
     endif else begin
        pinj[i]=0.
     endelse
  endfor
  
;=====================================
; Get injected Energy and particle mix
;=====================================
  power_mix1=fltarr(3)
  power_mix2=fltarr(3)
  energy1=fltarr(4)
  energy2=fltarr(4)
  physunit=0L
  error=0L
  dia_ref=0L
  time='                  '
  ed=0L
  exp='AUGD'
  type=2L
;; Open shotfile
  status=call_external(!libddww,'ddgetaug',   $
                       'ddopen',error,exp,diag,shot,ed,dia_ref,time)
  
;;Injector 1
  status=call_external(!libddww,'ddgetaug',   $
                       'ddparm',error,dia_ref,'INJ1    ','SPEC    ', $
                       type,3L,power_mix1,physunit)
  status=call_external(!libddww,'ddgetaug',   $
                       'ddparm',error,dia_ref,'INJ1    ','UEXQ    ', $
                       type,4L,energy1,physunit)
  
;;Injector 2
  status=call_external(!libddww,'ddgetaug',   $
                       'ddparm',error,dia_ref,'INJ2    ','SPEC    ', $
                       type,3L,power_mix2,physunit)
  status=call_external(!libddww,'ddgetaug',   $
                       'ddparm',error,dia_ref,'INJ2    ','UEXQ    ', $
                       type,4L,energy2,physunit)
  
;; close shotfile
  status=call_external(!libddww,'ddgetaug',   $
                       'ddclose',error,dia_ref)
  
;Injected Energy
  einj=[energy1,energy2]  
;Calculate Particle mix
  power_mix1=power_mix1*1.e-2
  power_mix2=power_mix2*1.e-2
  
   ;; print, ''
   ;; print, 'powermix box 1 (full, half, third):' $
   ;;        ,power_mix1(2), power_mix1(1), power_mix1(0)
   ;; print, 'powermix box 2 (full, half, third):' $
   ;;        , power_mix2(2), power_mix2(1), power_mix2(0)
  
  part_mix1=fltarr(3)
  part_mix2=fltarr(3)
  
  part_mix1(0)=power_mix1(2)/(power_mix1(2)+2.*power_mix1(1)+3.*power_mix1(0))
  part_mix1(1)=2.*power_mix1(1)/(power_mix1(2) $
                                 +2.*power_mix1(1) $
                                 +3.*power_mix1(0))
  part_mix1(2)=3.*power_mix1(0)/(power_mix1(2) $
                                 +2.*power_mix1(1) $
                                 +3.*power_mix1(0)) 
  
  part_mix2(0)=power_mix2(2)/(power_mix2(2) $
                              +2.*power_mix2(1)$
                              +3.*power_mix2(0))
  part_mix2(1)=2.*power_mix2(1)/(power_mix2(2) $
                                 +2.*power_mix2(1) $
                                 +3.*power_mix2(0))
  part_mix2(2)=3.*power_mix2(0)/(power_mix2(2) $
                                 +2.*power_mix2(1) $
                                 +3.*power_mix2(0)) 
  
  full=[part_mix1(0), part_mix1(0), part_mix1(0), part_mix1(0), $
        part_mix2(0), part_mix2(0), part_mix2(0), part_mix2(0) ]

  half=[part_mix1(1), part_mix1(1), part_mix1(1), part_mix1(1), $
        part_mix2(1), part_mix2(1), part_mix2(1), part_mix2(1) ]
  
  third=[part_mix1(2), part_mix1(2), part_mix1(2), part_mix1(2), $
         part_mix2(2), part_mix2(2), part_mix2(2), part_mix2(2) ]
 

  print, ''
  print, 'particle mix full:', full[[2,5]]
  print, 'particle mix half:', half[[2,5]]
  print, 'particle mix third:',third[[2,5]]
  print, ''
end
