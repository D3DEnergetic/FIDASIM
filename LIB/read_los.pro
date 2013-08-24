PRO read_los,file,los,save=save

	if file_test(file) then begin
       	ncid=ncdf_open(file,/nowrite)
		ncdf_varget,ncid,'Nchan',nchan
       	ncdf_varget,ncid,'xlens',xlens
       	ncdf_varget,ncid,'ylens',ylens
       	ncdf_varget,ncid,'zlens',zlens
       	ncdf_varget,ncid,'xlos',xlos
       	ncdf_varget,ncid,'ylos',ylos
       	ncdf_varget,ncid,'zlos',zlos
       	ncdf_varget,ncid,'headsize',headsize
       	ncdf_varget,ncid,'sigma_pi',sigma_pi
       	ncdf_varget,ncid,'opening_angle',opening_angle
       	ncdf_varget,ncid,'los_wght',weight

       	ncdf_close,ncid
		xyzlens=[[xlens],[ylens],[zlens]]
		xyzlos=[[xlos],[ylos],[zlos]]
		rlos=sqrt(xyzlos[*,0]^2+xyzlos[*,1]^2)
		
		los={nchan:nchan,rlos:rlos,xyzlens:xyzlens,xyzlos:xyzlos,$
			 headsize:headsize,opening_angle:opening_angle,sigma_pi:sigma_pi,weight:weight,err:0}
		if keyword_set(save) then save,los,filename='los.sav'
	endif else begin
		los={err:1}
	endelse
END
