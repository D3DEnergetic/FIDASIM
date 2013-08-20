PRO read_weights,file,weights,pitch_sign_convention=pitch_sign_convention
	
	if not keyword_set(pitch_sign_convention) then pitch_sign_convention=-1.0

	idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''

	if file_test(file) then begin
		
		;;READ IN ARRAY DIMENSION
		openr,55,file
	    readu,55 , idum & shot       = idum
	    readu,55 , fdum & time       = fdum 
	    readu,55 , idum & ichan_wght = idum
	    readu,55 , idum & nen        = idum
	    readu,55 , fdum & dE     = fdum
	    readu,55 , idum & npitch = idum
	    readu,55 , fdum & dpitch = fdum
	    readu,55 , idum & nwav   = idum
	    readu,55 , fdum & dwav   = fdum
	    readu,55 , fdum & wavel_start = fdum
	    readu,55 , fdum & wavel_end   = fdum
	    readu,55 , idum & nchan       = idum   

		;;DEFINE ARRAYS
	    angle=fltarr(nchan)
	    radius=fltarr(nchan)
	    weight_tot=fltarr(nchan,nwav,nen,npitch)
	    dummy_arr=fltarr(nwav,nen,npitch)

		;;READ IN ARRAYS
	    for chan=0,nchan-1 do begin
	       readu,55 , idum 
	       readu,55 , fdum & angle[chan]=fdum
    	   readu,55 , fdum & radius[chan]=fdum
	       readu,55, dummy_arr
	       weight_tot[chan,*,*,*]=dummy_arr
	    endfor
	    close,55

		;;ADDITIONALS PARAMETERS
	    central_wavel=(findgen(nwav)+0.5)*dwav+wavel_start
	    energyarr=(findgen(nen)+0.5)*dE
	    pitcharr=(findgen(npitch)+0.5)*dpitch-1.

	    pitcharr*=pitch_sign_convention
    	emax=max(energyarr)+0.5*dE

    	weights={nchan:nchan,ichan_wght:ichan_wght,nen:nen $
            ,dE:dE,emax:emax,emin:0.,npitch:npitch,dpitch:dpitch   $
            ,nwav:nwav,dwav:dwav,wavel_start:wavel_start,wavel_end:wavel_end   $
            ,central_wavel:central_wavel,energyarr:energyarr,pitcharr:pitcharr $
            ,weight_tot:weight_tot   $
            ,angle:angle,radius:radius,err:0}

	endif else begin
		weights={err:1}
	endelse

END
