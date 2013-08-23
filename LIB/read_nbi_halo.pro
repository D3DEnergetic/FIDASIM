PRO read_nbi_halo,file,nbi_halo,save=save

	if file_test(file) then begin
		ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'shot',shot
		ncdf_varget,ncid,'time',time
		ncdf_varget,ncid,'lambda',lambda
		ncdf_varget,ncid,'full',full
		ncdf_varget,ncid,'half',half
		ncdf_varget,ncid,'third',third
		ncdf_varget,ncid,'halo',halo
		ncdf_varget,ncid,'brems',brems
		ncdf_close,ncid

		nbi_halo={shot:shot,time:time,lambda:lambda,full:full,half:half,third:third,$
				  halo:halo,brems:brems,err:0}
		if keyword_set(save) then save,nbi_halo,filename='nbi_halo.sav'
	endif else begin
		nbi_halo={err:1}
	endelse

END
