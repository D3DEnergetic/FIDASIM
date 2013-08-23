PRO read_fida,file,fida,save=save

	if file_test(file) then begin

       	ncid=ncdf_open(file,/nowrite)
       	ncdf_varget,ncid,'shot',shot
       	ncdf_varget,ncid,'time',time
       	ncdf_varget,ncid,'lambda',lambda
       	ncdf_varget,ncid,'spectra',spectra
        ncdf_close,ncid

		fida={shot:shot,time:time,lambda:lambda,spectra:spectra,err:0}

		if keyword_set(save) then save,fida,filename='fida.sav'
	endif else begin
		fida={err:1}
	endelse
END
