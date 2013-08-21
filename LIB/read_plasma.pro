PRO read_plasma,plasma_file,plasma,save=save

	if file_test(plasma_file) then begin
		;;READ IN ARRAY SIZES
		nx=0L
		ny=0L
		nz=0L
		openr, lun, plasma_file, /get_lun
		readu,lun , nx
		readu,lun , ny
		readu,lun , nz

		;;DEFINE ARRAYS 
		bx=dblarr(nx,ny,nz)
		by=dblarr(nx,ny,nz)
		bz=dblarr(nx,ny,nz)
		ex=dblarr(nx,ny,nz)
		ey=dblarr(nx,ny,nz)
		ez=dblarr(nx,ny,nz)
		vti=dblarr(nx,ny,nz)
		dene=dblarr(nx,ny,nz)
		denp=dblarr(nx,ny,nz)
		denf=dblarr(nx,ny,nz)
		deni=dblarr(nx,ny,nz)
		te=dblarr(nx,ny,nz)
		ti=dblarr(nx,ny,nz)
		vrotx=dblarr(nx,ny,nz)
		vroty=dblarr(nx,ny,nz)
		vrotz=dblarr(nx,ny,nz)
		rho_grid=dblarr(nx,ny,nz)
		zeff=dblarr(nx,ny,nz)

		;;READ IN ARRAYS
		readu,lun , te
		readu,lun , ti
		readu,lun , dene
		readu,lun , denp
		readu,lun , deni
		readu,lun , denf
		readu,lun , vrotx
		readu,lun , vroty
		readu,lun , vrotz
		readu,lun , zeff
		readu,lun , bx
		readu,lun , by
		readu,lun , bz
		readu,lun , ex
		readu,lun , ey
		readu,lun , ez
		readu,lun , rho_grid
		close,lun
		free_lun, lun
		
		plasma={rho_grid:rho_grid,te:te,ti:ti,$
				dene:dene,deni:deni,denp:denp,denf:denf,zeff:zeff,$
				vrotx:vrotx,vroty:vroty,vrotz:vrotz,bx:bx,by:by,bz:bz,$
				ex:ex,ey:ey,ez:ez,err:0}
		if keyword_set(save) then save,plasma,filename='plasma.sav'
	endif else begin
		plasma={err:1}
	endelse

END
