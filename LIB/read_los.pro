PRO read_los,file,inputs,los

    idum=1L
    fdum=1.e0
    ddum=1.d0
    sdum=''

	if file_test(file) then begin

		openr, 55, file
		readu,55 , idum & nchan=idum
		xyzlos= FLTARR(nchan,3)
		xyzlens  = FLTARR(nchan,3)
		headsize = FLTARR(nchan) 
		opening_angle = FLTARR(nchan) 
		sigma_pi = FLTARR(nchan) 
		for i=0,nchan-1 do begin
			readu,55, ddum & xyzlens[i,0]=ddum
			readu,55, ddum & xyzlens[i,1]=ddum
			readu,55, ddum & xyzlens[i,2]=ddum
			readu,55, ddum & xyzlos[i,0]=ddum
			readu,55, ddum & xyzlos[i,1]=ddum
			readu,55, ddum & xyzlos[i,2]=ddum
			readu,55, ddum & headsize[i]=ddum
			readu,55, ddum & opening_angle[i]=ddum
			readu,55, ddum & sigma_pi[i]=ddum
		endfor
		weight=dblarr(inputs.nx,inputs.ny,inputs.nz,nchan)
		readu,55,weight
		close,55
	
		rlos=sqrt(xyzlos[*,0]^2+xyzlos[*,1]^2)
		
		los={nchan:nchan,rlos:rlos,xyzlens:xyzlens,xyzlos:xyzlos,$
			 headsize:headsize,opening_angle:opening_angle,sigma_pi:sigma_pi,weight:weight,err:0}
	endif else begin
		los={err:1}
	endelse
END
