@/afs/ipp/u/cxrs/idl/CFR_ausw/load_cfr.pro
pro cfr_setup,shot,det
  ;file='cfr_setup_'+string(shot,f='(i5)')+'.idl'
  ;if file_search(file) ne '' then begin
  ;   restore,file
  ;   return
  ;endif
  
  load_cfr,shot,cfr,/only_geom,err=err
  if err ne 0 then return
  nchan=cfr.header.ydim
  
  rlos=cfr.header.R_pos*100.d0  
  zlos=cfr.header.Z_pos*100.d0  
  philos=cfr.header.phi_pos 
  

  rhead=cfr.header.R_opt*100.d0  
  zhead=cfr.header.Z_opt*100.d0  
  phihead=cfr.header.phi_opt 

  ;Convert to uvz coordinates in [cm]
  ;At the intersection with NBI
  xlos=rlos*cos(philos)
  ylos=rlos*sin(philos)
  xhead = rhead*cos(phihead)
  yhead = rhead*sin(phihead)
   
  if shot gt 27500 then begin
    sigma_pi=replicate(0.9d0,nchan)
    index=where(zhead gt 100.,nind) 
    if nind gt 0 then sigma_pi[index]=1.  
  endif else begin
    sigma_pi=replicate(0.5d0,nchan)
  endelse
  det = set_det( NCHAN=nchan, $
                 XLOS=xlos, YLOS=ylos, ZLOS=zlos, $
                 XHEAD=xhead, YHEAD=yhead, ZHEAD=zhead, $
                 SIGMA_PI=sigma_pi )

  ;;save,filename=file,det

end
