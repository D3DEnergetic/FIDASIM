PRO read_birth,file,birth,save=save

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'birth_dens',birth_dens
		ncdf_close,ncid
	
		birth={shot:shot,time:time,birth_dens:birth_dens,err:0}
		
		if keyword_set(save) then save,birth,filename='birth.sav'
	endif else begin
		birth={err:1}
	endelse

END

