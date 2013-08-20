PRO read_fbm,file,FBM_struct,pitch_sign_convention=pitch_sign_convention
  ;---------------------------------------------------
  ; LOAD TRANSP fast ion distribution function
  ;---------------------------------------------------
  sdum=string(1,f='(i17)')
  ddum=1.d0
  idum=1L
  if not keyword_set(pitch_sign_convention) then pitch_sign_convention=-1.0

  ;;READ IN FILE
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
                    ,pitch_sign_convention:pitch_sign_convention,err:0}
  endif else begin
	FBM_struct={err:1}
  endelse

END

