FUNCTION newtfunc, xy 

common share, d,x210co,y210co,a210R,b210R

return, [(xy[1]-y210co)^2-d^2+(xy[0]-x210co)^2,$
	 xy[1]-b210R-a210R*xy[0]] 
 
END 

FUNCTION oblique_spatial,shot,whichcal=whichcal,plot=plot

;INPUT
;whichcal (optional)

;eventually will have new spatial calibrations
;for now,if whichcal is not specified, it uses last calibration from 2012

;OUTPUT
;returns x,y coords,major radius, and phi (usual DIII-D toroidal angle (deg)) of oblique FIDA spots at the midplane
;based on in-vessel spatial measurements

;CC 5/6/2014
;All 3 fibers per chord are theoretically supposed to lie along the
;same arc at the midplane
;So, the "radius" of each chord should be the mean value of the radius of
;the 3 fibers on the midplane.

;;;;modified based on /u/muscatel/fianalysis/2GFIDA/CALIB/spatial.pro


common share, d,x210co,y210co,a210R,b210R

IF not keyword_set(whichcal) THEN BEGIN
whichcal='forward2012'
print,'Using spatial calibration data from 4/10/2012'
ENDIF

;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
;for shots 
IF whichcal EQ 'forward2012' THEN BEGIN
caldate='20120410'
;;@@@@ 04/10/2012 calibrations @@@@
;;distance (cm) from crossover point to spot along target
;;cha   3      6      9      12    15    18    21    24    27     30     33;;;;
dspot1=[118.5, 112.5, 105.1, 99.0, 92.2, 86.0, 78.5, 72.0, 65.5, 59.2, 52.6]
;;cha   2      5      8      11    14    17    20    23    26     29     32;;;;
dspot2=[118.0, 112.5, 105.0, 98.5, 91.3, 85.2, 78.0, 71.0, 65.5, 57.6, 51.5]
;;cha   1      4      7      10    13    16    19    22    25     28     31;;;;
dspot3=[118.0, 112.0, 104.5, 98.0, 91.5, 84.6, 77.5, 70.8, 63.6, 57.0, 50.0]   
;dspot = (dspot1+dspot2+dspot3)/3.0
dspot = [dspot1,dspot2,dspot3]
;vertical distance (cm) from midplane to spot
zspot1=[8.0,    7.5,   6.0,   5.7,   5.3,   4.9,   4.3,   3.6,   3.0,   1.8,   1.3] 
zspot2=[0,     -1.4,  -2.0,  -3.4,  -4.0,  -5.6,  -7.2,  -7.5,  -8.8,  -9.2,  -10.0] 
zspot3=[-10.0, -11.2, -12.5, -15.0, -16.3, -17.0, -18.2, -19.0, -19.5, -20.5, -21.0] 
;zspot = (zspot1+zspot2+zspot3)/3.0
zspot = [zspot1,zspot2,zspot3]

fiber_index=[(indgen(11)*3+3),(indgen(11)*3+2),(indgen(11)*3+1)]
sorted=sort(fiber_index)

dspot=dspot(sorted)
zspot=zspot(sorted)

ENDIF

;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

nspot=n_elements(dspot)
xyrmid=fltarr(nspot,4) ;coords of spot [chord,x/y/R/phi] at midplane ;in cm

Rout=237.67 ;major radius (cm) of outer wall at midplane
Rin=101.40 ;major radius (cm) of inner wall at midplane

;from Scoville 2006 memo
t210Rsp=337.5 ;toroidal angle of 210R strikepoint
t210Rsp*=!pi/180. ;convert to radians

;From CER /u/antoniuk/b - stark/cerfit/source_cvs/cersrc/lib/rgnb/beam_geom.f90
;positive x is along the 90degree line, positive y is along the 0degree line 
x210co = -141.09; x coord (cm) of 210 crossover point
y210co = -230.40; y coord (cm) of 210 crossover point

;{x, y} coords of 210R strikepoint
x210Rsp = -Rout*Sin(2.*!pi - t210Rsp);
y210Rsp = Rout*Cos(2.*!pi - t210Rsp);

;slope of linear equation for 210R beamline
a210R = (y210Rsp - y210co)/(x210Rsp - x210co)

;y intercept of linear eq for beamlines
b210R = y210co - a210R*x210co;

;initial guess of X and Y points of spot on target
;  this is pretty forgiving and doesn't need to be changed
xy=[-127.,-111.] 

for i=0,nspot-1 do begin
  d=dspot[i]
  z=zspot[i]
  result=newton(xy,'newtfunc') ;{x,y} of spot on target (centerline of beam)
  xyztar=[result,z] ;{x,y,z} of spot through centerline of beam
  ;print,'{x,y,z} @ target=',xyztar
  
  ;**Now if z ne 0 then must find {x,y} of spot in midplane***

  detp=100.*[-0.4602,-1.9850,1.2200] ;{x,y,z} of 195V+3 optical assembly

  ;solve for {x,y,R} at midplane (z=0)
  zmid=0.
  xyrmid[i,0]=(zmid-xyztar[2])*(detp[0]-xyztar[0])/(detp[2]-xyztar[2])+xyztar[0] ;x
  xyrmid[i,1]=(zmid-xyztar[2])*(detp[1]-xyztar[1])/(detp[2]-xyztar[2])+xyztar[1] ;y
  xyrmid[i,2]=sqrt(xyrmid[i,0]^2+xyrmid[i,1]^2) ;R
  xyrmid[i,3]=abs(atan(xyrmid[i,0]/xyrmid[i,1]))*180./!pi+180. ;phi, usual DIII-D toroidal angle (deg)


end ;for

;get mean value R for each chord
meanRmid=fltarr(11) ;mean R at midplane;in cm
meanPhimid=fltarr(11) ;mean Phi at midplane
xx=fltarr(11)
yy=fltarr(11)

for i=0,10 do begin
meanRmid[i]= mean([xyrmid[i*3,2],xyrmid[i*3+1,2],xyrmid[i*3+2,2]]) ;mean R
meanPhimid[i]= mean([xyrmid[i*3,3],xyrmid[i*3+1,3],xyrmid[i*3+2,3]]) ;mean Phi
xx[i]=meanRmid[i]*sin(meanPhimid[i]*!pi/180)
yy[i]=meanRmid[i]*cos(meanPhimid[i]*!pi/180)

end ;for




if keyword_set(plot) then begin
  loadct,39,/SILENT

  !p.charsize=1.5 & !p.charthick=5
  !p.font=0
  xout=(2.*findgen(100)/99.-1.)*rout
  xin=(2.*findgen(100)/99.-1.)*rin
  plot,xout,sqrt(rout^2-xout^2),xrange=[-rout,0],yrange=[-rout,0],$
    thick=4,/iso,xtitle='X (cm)',ytitle='Y (cm)',color=1, background=255
  oplot,xout,-sqrt(rout^2-xout^2),thick=4,color=1
  oplot,xin,sqrt(rin^2-xin^2),thick=4,color=1
  oplot,xin,-sqrt(rin^2-xin^2),thick=4,color=1
  xbeam=findgen(100)/99.*(-100.)-100
  oplot,xbeam,a210R*xbeam+b210R,thick=12,color=50
  oplot,xyrmid[*,0],xyrmid[*,1],psym=4,color=240,symsize=1,thick=3
  

oplot,xx,yy,psym=5,color=100,symsize=1,thick=2

end ;plot

spatial_data=CREATE_STRUCT('caldate',caldate)

FOR i=0,10 DO BEGIN

fibers_xyr={fibers:[i*3+1,i*3+2,i*3+3],$
            x:[xyrmid[i*3,0],xyrmid[i*3+1,0],xyrmid[i*3+2,0]],$
            y:[xyrmid[i*3,1],xyrmid[i*3+1,1],xyrmid[i*3+2,1]],$
            r:[xyrmid[i*3,2],xyrmid[i*3+1,2],xyrmid[i*3+2,2]],$
          phi:[xyrmid[i*3,3],xyrmid[i*3+1,3],xyrmid[i*3+2,3]]}

chord_struc = {X:xx[i],$
               Y:yy[i],$
               R:meanRmid[i],$
               Phi:meanPhimid[i],$
               Fibers:fibers_xyr,$
               Units:'cm and degrees'}

spatial_data=CREATE_STRUCT(spatial_data,'P'+STRTRIM(i+1,2),chord_struc)


ENDFOR

return,spatial_data

END

