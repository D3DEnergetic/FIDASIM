PRO read_nbi_halo,file,nbi_halo,save=save

    idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''

	if file_test(file) then begin
		;;READ IN ARRAY SIZES
		openr,55,file
		readu,55,idum  ;; shot
		readu,55,fdum  ;; time
		readu,55,idum & nlos=idum 
		readu,55,idum & nlambda=idum

		;;DEFINE ARRAYS
		lambda=fltarr(nlambda)
		fspectra=fltarr(nlambda,nlos)
		hspectra=fltarr(nlambda,nlos)
		tspectra=fltarr(nlambda,nlos)
		halospectra=fltarr(nlambda,nlos)
		bremspectra=fltarr(nlambda,nlos)

		;;READ IN ARRAYS
		readu,55,lambda
		readu,55,fspectra
		readu,55,hspectra
		readu,55,tspectra
		readu,55,halospectra
		readu,55,bremspectra
		close,55
		
		nbi_halo={lambda:lambda,full:fspectra,half:hspectra,third:tspectra,$
				  halo:halospectra,brems:bremspectra,err:0}
		if keyword_set(save) then save,nbi_halo,filename='nbi_halo.sav'
	endif else begin
		nbi_halo={err:1}
	endelse

END
