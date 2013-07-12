PRO d3d_plots,inputs,grid, nbi, fida, equil,rot_mat

g=equil.g

	;;PLOTTING PLANE VIEW BEAMS AND CHORDS
	window,0 & wset,0
	loadct,39,/silent

	;;GET PROPER RANGES
	xmin=min(grid.xx) & ymin=min(grid.yy) & zmin=min(grid.zz)
	xmax=max(grid.xx) & ymax=max(grid.yy) & zmax=max(grid.zz)
	if xmin lt 0 then xmin1=1.2*xmin else xmin1=0.8*xmin
	if xmax lt 0 then xmax1=0.8*xmax else xmax1=1.2*xmax
	if ymin lt 0 then ymin1=1.2*ymin else ymin1=0.8*ymin
	if ymax lt 0 then ymax1=0.8*ymax else ymax1=1.2*ymax
	if zmin lt 0 then zmin1=1.2*zmin else zmin1=0.8*zmin
	if zmax lt 0 then zmax1=0.8*zmax else zmax1=1.2*zmax
	x_range=[xmin1,xmax1] & y_range=[ymin1,ymax1] & z_range=[zmin1,zmax1]

x_range=[xmin1,xmax1] & y_range=[ymin1,ymax1]

	plot,[0],[0],/nodata,xrange=x_range,yrange=y_range,$
		color=0,background=255,title='PLAN VIEW',xtitle='X [cm]',ytitle='Y [cm]'
	for i=0,grid.nx-1 do begin
		oplot,[grid.xx[i],grid.xx[i]],[grid.yy[0],grid.yy[-1]],color=0
	endfor	
	for i=0,grid.ny-1 do begin
		oplot,[grid.xx[0],grid.xx[-1]],[grid.yy[i],grid.yy[i]],color=0
	endfor	

	for i=0,fida.nchan-1 do begin
		oplot,[fida.xmid[i],fida.xlens[i]] ,[fida.ymid[i],fida.ylens[i]] ,color=50
	endfor
	ii=inputs.isource[0]
	uvw_ray=[-1.d0,0.d0,0.d0]*1000.
    updown=1
    rotate_uvw,uvw_ray,rot_mat.Arot,rot_mat.Brot,rot_mat.Crot,updown $
              ,xyz_ray
    oplot,[nbi.xyz_src[ii,0],nbi.xyz_src[ii,0]+xyz_ray[0]] $
         ,[nbi.xyz_src[ii,1],nbi.xyz_src[ii,1]+xyz_ray[1]],thick=2,color=230

w=where(g.bdry[0,*] gt 0.)
rmin=100.*min(g.bdry[0,w]) & rmax=100.*max(g.bdry[0,w])
rmaxis=100.*g.rmaxis
phi=2.*!pi*findgen(501)/500.
oplot,rmin*cos(phi),rmin*sin(phi),color=150
oplot,rmaxis*cos(phi),rmaxis*sin(phi),color=150,linestyle=2
oplot,rmax*cos(phi),rmax*sin(phi),color=150

;----------------------------------------------
	;;PLOT CROSS SECTION BEAM AND CHORDS 
	window,1 & wset,1
	plot,[0],[0],/nodata,xrange=[rmin,rmax], $
            yrange=100.*[min(g.lim[1,*]),max(g.lim[1,*])],$
		color=0,background=255,title='ELEVATION',xtitle='R [cm]',ytitle='Z [cm]'

oplot,grid.r_grid,grid.zc,psym=3  
; Lines of sight
for i=0,fida.nchan-1 do begin
  if fida.zlens[i] ne fida.zmid[i] then begin
    z=(fida.zlens[i]-fida.zmid[i])*findgen(101)/100.+fida.zmid[i]
    x=(fida.xlens[i]-fida.xmid[i])*(z-fida.zmid[i])/ $
      (fida.zlens[i]-fida.zmid[i]) + fida.xmid[i]
    y=(fida.ylens[i]-fida.ymid[i])*(z-fida.zmid[i])/ $
      (fida.zlens[i]-fida.zmid[i]) + fida.ymid[i]
    oplot,sqrt(x^2+y^2),z,color=50
  end else begin 
    y=(fida.ylens[i]-fida.ymid[i])*findgen(101)/100.+fida.ymid[i]
    x=(fida.xlens[i]-fida.xmid[i])*(y-fida.ymid[i])/ $
      (fida.ylens[i]-fida.ymid[i]) + fida.ymid[i]
    oplot,sqrt(x^2+y^2),replicate(fida.zlens[i],101),color=50
  end
end
; Equilibrium
oplot,100.*g.bdry[0,*],100.*g.bdry[1,*],color=150
oplot,100.*g.lim[0,*],100.*g.lim[1,*],color=0

end
