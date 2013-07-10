@/afs/ipp/u/cxrs/idl/CFR_ausw/load_cfr.pro
pro cfr_setup,shot,det
  
  load_cfr,shot,cfr,/only_geom
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
  det={ nchan:nchan $
        ,xlos:xlos, ylos:ylos, zlos: zlos $
        ,xhead: xhead, yhead: yhead, zhead: zhead, headsize: fltarr(nchan)}

  save,filename=file,det

end
