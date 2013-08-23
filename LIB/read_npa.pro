PRO read_npa,file,npa,save=save

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'ipos',npaipos
		ncdf_varget,ncid,'fpos',npafpos
		ncdf_varget,ncid,'v',npav
		ncdf_varget,ncid,'wght',npawght
		ncdf_close,ncid

		if n_elements(npawght) gt 0 then begin
			npa={npaipos:npaipos,npafpos:npafpos,npav:npav, npawght:npawght,err:0}
			if keyword_set(save) then save,npa,filename='npa.sav'
		endif else begin
			print, 'attention, no particle reached the detector!'
			npa={err:1}
		endelse
	endif else begin
		npa={err:1}
	endelse
END
