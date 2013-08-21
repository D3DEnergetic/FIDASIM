PRO read_npa,file,npa,save=save

	if file_test(file) then begin

		openr, 55, npa_file
		readu,55,idum;& print, 'shot:', fdum
		readu,55,fdum;& print, 'time:', fdum
		readu,55,idum & nr_npa=idum
		readu,55,idum & counter=idum

		if counter gt 0 then begin
			npaipos=fltarr(counter,3)
			npafpos=fltarr(counter,3)
			npav  = fltarr(counter,3)
			npawght=fltarr(counter)
			readu,55,npaipos
			readu,55,npafpos
			readu,55,npav
			readu,55,npawght 
			close,55
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
