PRO read_birth,file,birth,save=save

    idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''

	if file_test(file) then begin
		;;READ IN ARRAY SIZES
		openr,55,file
		readu,55,fdum & shot=fdum
		readu,55,fdum & time=fdum
		readu,55,fdum & nx=fdum
		readu,55,fdum & ny=fdum
		readu,55,fdum & nz=fdum
		readu,55,fdum & npitch=fdum

		;;DEFINE ARRAYS
		birth_dens=fltarr(nx,ny,nz,3,npitch)
		
		;;READ IN ARRAYS
		readu,55,birth_dens
		close,55
	
		birth={birth_dens:birth_dens,err:0}
		
		if keyword_set(save) then save,birth,filename='birth.sav'
	endif else begin
		birth={err:1}
	endelse

END

