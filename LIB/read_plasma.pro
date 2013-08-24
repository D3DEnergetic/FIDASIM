PRO read_plasma,plasma_file,plasma,save=save

	if file_test(plasma_file) then begin

       	ncid=ncdf_open(plasma_file,/nowrite)
       	ncdf_varget,ncid,'ti',ti
       	ncdf_varget,ncid,'te',te
       	ncdf_varget,ncid,'dene',dene
       	ncdf_varget,ncid,'deni',deni
       	ncdf_varget,ncid,'denp',denp
       	ncdf_varget,ncid,'denf',denf
       	ncdf_varget,ncid,'vrotx',vrotx
       	ncdf_varget,ncid,'vroty',vroty
       	ncdf_varget,ncid,'vrotz',vrotz
       	ncdf_varget,ncid,'zeff',zeff
       	ncdf_varget,ncid,'bx',bx
       	ncdf_varget,ncid,'by',by
       	ncdf_varget,ncid,'bz',bz
       	ncdf_varget,ncid,'ex',ex
       	ncdf_varget,ncid,'ey',ey
       	ncdf_varget,ncid,'ez',ez
       	ncdf_varget,ncid,'rho_grid',rho_grid

       	ncdf_close,ncid	
		
		plasma={rho_grid:rho_grid,te:te,ti:ti,$
				dene:dene,deni:deni,denp:denp,denf:denf,zeff:zeff,$
				vrotx:vrotx,vroty:vroty,vrotz:vrotz,bx:bx,by:by,bz:bz,$
				ex:ex,ey:ey,ez:ez,err:0}
		if keyword_set(save) then save,plasma,filename='plasma.sav'
	endif else begin
		plasma={err:1}
	endelse

END
