PRO read_neutrals,file,neutrals,save=save

    idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''

	if file_test(file) then begin
		;;READ IN ARRAY SIZES
    	openr, 55, file
     	readu,55,dum
     	readu,55,dum
     	readu,55,dum  & nx        = fix(dum)  
     	readu,55,dum  & ny        = fix(dum)   
     	readu,55,dum  & nz        = fix(dum)    
     	readu,55,dum  & nlevs     = fix(dum) 

		;;DEFINE ARRAYS
     	fdens   =fltarr(nx,ny,nz,nlevs)
     	hdens   =fltarr(nx,ny,nz,nlevs)
     	tdens   =fltarr(nx,ny,nz,nlevs)
     	halodens=fltarr(nx,ny,nz,nlevs)
     	
		;;READ IN ARRAYS
     	readu,55,fdens
     	readu,55,hdens
     	readu,55,tdens
     	readu,55,halodens
     	close,55

     	neutrals={fdens:fdens, hdens:hdens,tdens:tdens, halodens:halodens,err:0}
		if keyword_set(save) then save,neutrals,filename='neutrals.sav'
	endif else begin
		neutrals={err:1}
  	endelse

END
