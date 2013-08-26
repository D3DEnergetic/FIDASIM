PRO augd_routines,inputs,grid,nbi,chords,profiles,equil,err
  user=+GETENV('USER')

  ;;GET PROFILES
;  load_transp_profiles,inputs.shot,inputs.time,profiles,inputs.cdf_file
  restore,'/u/'+user+'/FIDASIM/TEST/load_transp_profiles.idl'
  profiles=iprof

  ;;GET B FIELD
;  load_bfield,inputs,grid,bx,by,bz,rhopf,rhotf
;  if profiles.rho_str eq 'rho_tor' then rho_grid=rhotf 
;  if profiles.rho_str eq 'rho_pol' then rho_grid=rhopf  
  restore,'/u/'+user+'/FIDASIM/TEST/load_bfield.idl'

  equil={rho_grid:rho_grid,rho_chords:{rhos:0,ds:0},bx:bx,by:by,bz:bz,ex:bx*0,ey:by*0,ez:bz*0}

  ;; NBI GEOMETRY and DATA
  ;nbi_data,inputs,einj,pinj $ ;power,energy,particle mix
  ;              ,ffull,fhalf,fthird,doplot=0 $
  ;              ,ps=inputs.ps  
  restore,'/u/'+user+'/FIDASIM/TEST/nbi_data.idl'
  nbi_geometry_transp,nbgeom
  help,nbgeom
  nbi={einj:double(einj[inputs.isource]),pinj:pinj[inputs.isource],full:double(ffull[inputs.isource]),$
	   half:double(fhalf[inputs.isource]),third:double(fthird[inputs.isource]),$
	   xyz_src:reform(nbgeom.xyz_src[inputs.isource,*]),xyz_pos:reform(nbgeom.xyz_pos[inputs.isource,*]),$
	   bmwidra:nbgeom.bmwidra,bmwidza:nbgeom.bmwidza,$
	   divy:nbgeom.divy[*,inputs.isource],divz:nbgeom.divz[*,inputs.isource],focy:nbgeom.focy[inputs.isource],focz:nbgeom.focz[inputs.isource]}

  ;;GET CHORD GEOMETRY
  if strupcase(inputs.diag) eq 'NPA' then begin
      npa_setup,det
  endif else begin
;	cfr_setup,inputs.shot,det
	restore,'/u/'+user+'/FIDASIM/TEST/cfr_setup.idl'
  endelse

  chords={sigma_pi_ratio:det.sigma_pi,nchan:det.nchan,xlos:det.xlos,ylos:det.ylos,$
		zlos:det.zlos,xlens:det.xhead,ylens:det.yhead,zlens:det.zhead,$
		opening_angle:det.npa.opening_angle,headsize:det.npa.headsize}

 err=0
END
