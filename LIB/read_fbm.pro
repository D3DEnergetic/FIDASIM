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

