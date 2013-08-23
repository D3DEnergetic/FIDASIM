PRO read_neutrals,file,neutrals,save=save

	if file_test(file) then begin
		ncid=NCDF_OPEN(file,/nowrite)
		
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'fdens',fdens
		ncdf_varget,ncid,'hdens',hdens
		ncdf_varget,ncid,'tdens',tdens
		ncdf_varget,ncid,'halodens',halodens
		NCDF_CLOSE,ncid

     	neutrals={shot:shot,time:time,fdens:fdens, hdens:hdens,tdens:tdens, halodens:halodens,err:0}
		if keyword_set(save) then save,neutrals,filename='neutrals.sav'
	endif else begin
		neutrals={err:1}
  	endelse

END
