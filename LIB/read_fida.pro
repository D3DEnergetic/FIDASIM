PRO read_fida,file,fida,save=save

	idum=1L
	fdum=1.e0
	ddum=1.d0
	sdum=''

	if file_test(file) then begin
		;;READ IN ARRAY SIZES
		openr,55,file
		readu,55,idum ;;SHOT
		readu,55,fdum ;;TIME
		readu,55,idum & nlos=idum
		readu,55,idum & nlambda=idum
		
		;;DEFINE ARRAYS
		lambda=fltarr(nlambda)
		spectra=fltarr(nlambda,nlos)

		;;READ IN ARRAYS
		readu,55,lambda
		readu,55,spectra
		close,55

		fida={lambda:lambda,spectra:spectra,err:0}

		if keyword_set(save) then save,fida,filename='fida.sav'
	endif else begin
		fida={err:1}
	endelse

END
