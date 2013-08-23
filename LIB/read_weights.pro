PRO read_weights,file,weights,save=save,pitch_sign_convention=pitch_sign_convention
	
	if not keyword_set(pitch_sign_convention) then pitch_sign_convention=-1.0

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'lambda',central_wavel
		ncdf_varget,ncid,'energy',energyarr
		ncdf_varget,ncid,'pitch',pitcharr
		ncdf_varget,ncid,'radius',radarr
		ncdf_varget,ncid,'theta',theta_arr
		ncdf_varget,ncid,'wfunct',weight_tot
		ncdf_close,ncid
		
		dE=abs(energyarr[1]-energyarr[0])
		dpitch=abs(pitcharr[1]-pitcharr[0])
		dwav=abs(central_wavel[1]-central_wavel[0])

		nwav=n_elements(central_wavel)
	    nen=n_elements(energyarr)
		npitch=n_elements(pitcharr)

	    pitcharr*=pitch_sign_convention
    	emax=max(energyarr)+0.5*dE

    	weights={shot:shot,time:time,nchan:n_elements(radarr),nen:nen $
            ,dE:dE,emax:emax,emin:0.,npitch:npitch,dpitch:dpitch   $
            ,nwav:nwav,dwav:dwav   $
            ,central_wavel:central_wavel,energyarr:energyarr,pitcharr:pitcharr $
            ,weight_tot:weight_tot   $
            ,angle:theta_arr,radius:radarr,err:0}
		if keyword_set(save) then save,weights,filename='weights.sav'
	endif else begin
		weights={err:1}
	endelse

END
