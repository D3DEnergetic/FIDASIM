PRO visualize_fidasim_inputs,file_spec,noplasma=noplasma
; Plot FIDASIM input files
; WWH 9/14/2016

; INPUT
; file_spec  prefix of input files, e.g.,
;    '/p/fida/FIDASIM/RESULTS/D3D/153072H01_2016/3300_'

; OUTPUT
;  elevation contour plots of density, Te, Ti, v_toroidal, Zeff, denf
;  plan & elevation of grid, beam, and sightlines

; KEYWORD
; noplasma  skip plots of plasma parameters and beam density

; USES
; .compile /p/fida/FIDASIM/lib/read_hdf5
; .compile /p/fida/FIDASIM/lib/readg

; ASSUMES: 
;  an EFIT equilibrium
;  on-axis active beam

;--------------------
; Check that files exist
zz=findfile(file_spec+'*')
if zz[0] eq '' then print, $
  'NO FILES FOUND:',file_spec

d=read_hdf5(file_spec+'distribution.h5',/shallow)
e=read_hdf5(file_spec+'equilibrium.h5',/shallow)
geo=read_hdf5(file_spec+'geometry.h5',/shallow)
nbi=geo.nbi
if nbi.axis ne 0. or nbi.aoffz ne 0. then print, $
  '*** WARNING! Program does not plot off-axis beams correctly', $
  nbi.axis,nbi.aoffz
g=readg(e.fields.data_source)
if g.error eq 1 then print, $
  'EQDSK NOT FOUND:',e.fields.data_source

;---------------
; Extract grid information from inputs.dat 
openr,inunit,file_spec+'inputs.dat', /get_lun
line=' '
while(not eof(inunit)) do begin
      readf,inunit,line
      CASE 1 of
      strpos(line,'xmin') gt -1: begin
        xp=strpos(line,'!!')
        xmin=float(strmid(line,7,xp-8))
      end
      strpos(line,'xmax') gt -1: begin
        xp=strpos(line,'!!')
        xmax=float(strmid(line,7,xp-8))
      end
      strpos(line,'ymin') gt -1: begin
        xp=strpos(line,'!!')
        ymin=float(strmid(line,7,xp-8))
      end
      strpos(line,'ymax') gt -1: begin
        xp=strpos(line,'!!')
        ymax=float(strmid(line,7,xp-8))
      end
      strpos(line,'zmin') gt -1: begin
        xp=strpos(line,'!!')
        zmin=float(strmid(line,7,xp-8))
      end
      strpos(line,'zmax') gt -1: begin
        xp=strpos(line,'!!')
        zmax=float(strmid(line,7,xp-8))
      end
      strpos(line,'alpha') gt -1: begin
        xp=strpos(line,'!!')
        alpha=float(strmid(line,8,xp-9))
      end
      strpos(line,'beta') gt -1: begin
        xp=strpos(line,'!!')
        beta=float(strmid(line,8,xp-9))
      end
      strpos(line,'gamma') gt -1: begin
        xp=strpos(line,'!!')
        gamma=float(strmid(line,8,xp-9))
      end
      strpos(line,'origin(1)') gt -1: begin
        origin=fltarr(3)
        xp=strpos(line,'!!')
        origin[0]=float(strmid(line,12,xp-13))
      end
      strpos(line,'origin(2)') gt -1: begin
        xp=strpos(line,'!!')
        origin[1]=float(strmid(line,12,xp-13))
      end
      strpos(line,'origin(3)') gt -1: begin
        xp=strpos(line,'!!')
        origin[2]=float(strmid(line,12,xp-13))
      end
      strpos(line,'calc_npa_wght') gt -1: begin
        xp=strpos(line,'!!')
        calc_npa_wght=fix(strmid(line,16,xp-17))
      end
      strpos(line,'calc_fida_wght') gt -1: begin
        xp=strpos(line,'!!')
        calc_fida_wght=fix(strmid(line,17,xp-18))
      end
      strpos(line,'calc_fida') gt -1: begin
        xp=strpos(line,'!!')
        calc_fida=fix(strmid(line,12,xp-13))
      end
      strpos(line,'calc_npa') gt -1: begin
        xp=strpos(line,'!!')
        calc_npa=fix(strmid(line,11,xp-12))
      end
      else:
  endcase
  endwhile
close,/all
help,xmin,xmax,ymin,ymax,zmin,zmax,alpha,beta,gamma
print,origin
if beta ne 0. or gamma ne 0. then print, $
 'WARNING!  beta and/or gamma are non-zero:',beta,gamma
uorigin=origin[0]+xmin*cos(alpha)
vorigin=origin[1]+xmin*sin(alpha)
umax=origin[0]+xmax*cos(alpha)
vmax=origin[1]+xmax*sin(alpha)
rorigin=sqrt(uorigin^2+vorigin^2)
rmax=sqrt(umax^2+vmax^2)

dospec=calc_fida>calc_fida_wght
donpa=calc_npa>calc_npa_wght
if dospec eq 0 then print, $
   'FIDA will NOT be calculated' else $
  spec=geo.spec

if donpa eq 0 then print, $
   'NPA will NOT be calculated' else $
  npa=geo.npa

;--------------------
; Plasma parameters

loadct,13
if ~keyword_set(noplasma) then begin

window,/free
a=e.plasma.te
contour,a,e.plasma.r,e.plasma.z,nlevels=124,/fill, $
  title='Te '+string(min(a))+string(max(a))+' keV', $
  xtitle='R (cm)',ytitle='z (cm)'
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2

window,/free
a=e.plasma.ti
contour,a,e.plasma.r,e.plasma.z,nlevels=124,/fill, $
  title='Ti '+string(min(a))+string(max(a))+' keV', $
  xtitle='R (cm)',ytitle='z (cm)'
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2

window,/free
a=e.plasma.dene
contour,a,e.plasma.r,e.plasma.z,nlevels=124,/fill, $
  title='Density '+string(min(a))+string(max(a))+' cm^-3', $
  xtitle='R (cm)',ytitle='z (cm)'
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2

window,/free
a=e.plasma.vt
contour,a,e.plasma.r,e.plasma.z,nlevels=124,/fill, $
  title='Toroidal rotation '+string(min(a))+string(max(a))+' cm/s', $
  xtitle='R (cm)',ytitle='z (cm)'
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2

window,/free
a=e.plasma.zeff
contour,a,e.plasma.r,e.plasma.z,nlevels=124,/fill, $
  title='Zeff '+string(min(a))+string(max(a)), $
  xtitle='R (cm)',ytitle='z (cm)'
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2

window,/free
a=d.denf
contour,a,e.plasma.r,e.plasma.z,nlevels=124,/fill, $
  title='Beam density '+string(min(a))+string(max(a))+' cm^-3', $
  xtitle='R (cm)',ytitle='z (cm)'
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2

end  ; doplasma

;-----------------
; Sightlines
; Plan view

; Plot grid first
window,/free,xsize=600,ysize=600
bdrymax=max(g.bdry[0,*])*100
bdrymin=min(g.bdry[0,0:g.nbdry-1])*100
bdryzmax=max(g.bdry[1,*])*100
bdryzmin=min(g.bdry[1,*])*100

rplot=bdrymax>rorigin
plot,[0],xrange=[-rplot,rplot],yrange=[-rplot,rplot], $
  /xstyle,/ystyle,/nodata
theta=2*!pi*findgen(501)/500.
oplot,bdrymax*cos(theta),bdrymax*sin(theta)
oplot,bdrymin*cos(theta),bdrymin*sin(theta)
oplot,[uorigin,umax],[vorigin,vmax],thick=3,color=200
oplot,[uorigin,umax]-ymax*sin(alpha),[vorigin,vmax]+ymax*cos(alpha), $
  linestyle=2,color=200
oplot,[uorigin,umax]-ymin*sin(alpha),[vorigin,vmax]+ymin*cos(alpha), $
  linestyle=2,color=200

; Active Beam
xnbi=nbi.src[0] + nbi.axis[0]*findgen(2000)
ynbi=nbi.src[1] + nbi.axis[1]*findgen(2000)
rnbi=sqrt(xnbi^2+ynbi^2)
wnbi=where(rnbi ge bdrymin)
xnbi=xnbi[wnbi] & ynbi=ynbi[wnbi] & rnbi=rnbi[wnbi]
oplot,xnbi,ynbi,color=150
anbi=atan(nbi.axis[1],nbi.axis[0])
du=nbi.awidy[0]*sin(anbi)
dv=nbi.awidy[0]*cos(anbi)
oplot,xnbi+du,ynbi+dv,linestyle=2,color=150
oplot,xnbi-du,ynbi-dv,linestyle=2,color=150

s=300.*findgen(101)/100.
uvz=fltarr(3,101)
; FIDA
if dospec then begin
 for j=0,spec.nchan-1 do begin
 for i=0,2 do uvz[i,*]=spec.lens[i,j]+s*spec.axis[i,j]
 rr=sqrt(uvz[0,*]^2+uvz[1,*]^2)
 w=where(rr ge bdrymin and rr le bdrymax and $
        uvz[2,*] ge bdryzmin and uvz[2,*] le bdryzmax,nw) 
 if nw gt 0 then oplot,uvz[0,w],uvz[1,w],color=100 else $
   print,'Chord outside box:',j,spec.radius[j]
 end
end

; NPA
if donpa then begin
  axis=fltarr(3)
  for j=0,npa.nchan-1 do begin
    axis=reform(npa.a_cent[*,j]-npa.d_cent[*,j])
    axis/=sqrt(axis[0]^2+axis[1]^2+axis[2]^2)
    for i=0,2 do uvz[i,*]=npa.d_cent[i,j]+s*axis[i]
    rr=sqrt(uvz[0,*]^2+uvz[1,*]^2)
    w=where(rr ge bdrymin and rr le bdrymax and $
        uvz[2,*] ge bdryzmin and uvz[2,*] le bdryzmax,nw) 
    if nw gt 0 then oplot,uvz[0,w],uvz[1,w],color=60 else $
      print,'Chord outside box:',j,spec.radius[j]
  end
end

; Elevation
zpmin=bdryzmin<zmin
zpmax=bdryzmax>zmax
rpmin=min(100*g.lim[0,*])<rmax
rpmax=max(100*g.lim[0,*])>rorigin
window,/free,xsize=500,ysize=fix(500*(zpmax-zpmin)/(rpmax-rpmin))

plot,[0],/nodata,xrange=[rpmin,rpmax],yrange=[zpmin,zpmax],/xstyle,/ystyle
oplot,100*g.bdry[0,*],100*g.bdry[1,*],thick=2
oplot,100*g.lim[0,*],100*g.lim[1,*],thick=2
oplot,[rorigin,rmax],replicate(0.5*(zmin+zmax),2),thick=3,color=200
oplot,[rorigin,rmax],replicate(zmin,2),linestyle=2,color=200
oplot,[rorigin,rmax],replicate(zmax,2),linestyle=2,color=200
oplot,[min(rnbi),max(rnbi)],replicate(nbi.aoffz,2),color=150
oplot,[min(rnbi),max(rnbi)], $
          replicate(nbi.aoffz+nbi.awidz,2),linestyle=2,color=150
oplot,[min(rnbi),max(rnbi)], $
          replicate(nbi.aoffz-nbi.awidz,2),linestyle=2,color=150

if dospec then begin
 for j=0,spec.nchan-1 do begin
  for i=0,2 do uvz[i,*]=spec.lens[i,j]+s*spec.axis[i,j]
  rr=sqrt(uvz[0,*]^2+uvz[1,*]^2)
  w=where(rr ge bdrymin and rr le bdrymax and $
        uvz[2,*] ge bdryzmin and uvz[2,*] le bdryzmax,nw) 
  if nw gt 0 then oplot,rr[w],uvz[2,w],color=100
 end
end

if donpa then begin
  axis=fltarr(3)
  for j=0,npa.nchan-1 do begin
    axis=reform(npa.a_cent[*,j]-npa.d_cent[*,j])
    axis/=sqrt(axis[0]^2+axis[1]^2+axis[2]^2)
    for i=0,2 do uvz[i,*]=npa.d_cent[i,j]+s*axis[i]
    rr=sqrt(uvz[0,*]^2+uvz[1,*]^2)
    w=where(rr ge bdrymin and rr le bdrymax and $
        uvz[2,*] ge bdryzmin and uvz[2,*] le bdryzmax,nw) 
    if nw gt 0 then oplot,rr[w],uvz[w,w],color=60
  end
end


end
