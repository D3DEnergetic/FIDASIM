;***********************************************************************************
;
;  4dlib/EFITLIB/kupfer/idl_geqdsk/contour_psi.pro
;
;  created:
;
;  modified:
;   7/25/01     Ted Terpstra - put in  psi =  temporary(psi[i])
;               to fix a problem for Degrassi who had a function psi.pro.
;   3/02/99     Qian Peng - updates for efit
;   8/06/98	Jeff Schachter - removed call to CONVERT_COORD
;		(unnecessary, and as of IDL 5.1, illegal)
;  10/08/97     MRW - To avoid all the problem with the path.dat file,
;                  changed algorithm to not use file output - use
;                  path_info and path_xy in call to contour instead
;                  of path_filename.
;   9/23/97     TBT - write path.dat file to /scratch as per TOdd Evans
;   8/12/97     TBT - put in print of filename to file path.dat
;                  There was a conflict with the path.dat filename
;                  because two jobs were calling this at the same time.
;                  Added hr, mn sc to filename to make it
;                     pathhrmnsc.dat where 00<=hr<24  0<=mn<60 0<=sc<60
;  10/31/96     Ted Terpstra - Happy Halloween
;               Solved Dan Baker's problem by setting scaling
;               back to default for call to contour.
;    8-9-96     K. Greene -- Add use of function TEMPORARY
;    6-7-96	K. Greene -- Reformat source file
;
;
;***********************************************************************************
;
; INPUTS -- 
;     	g = data structure containing GEQDSK (from readg routine).
;     	npts = number of theta points
;     	psi1 = normalized value of first psi surface.
;     	psi2 = normalized value of last psi surface.
;     	nx = the number of x points on mesh (minimum value nx=g.mw)
;     	ny = number of y points on mesh (minimum value ny=g.mh)
; 		Here "mesh" refers to the cartesian (x,y) mesh where psi is
; 		tabulated. The minimum values of nx and ny correspond to 
; 		choosing the original mesh specified in the GEQDSK.
; 		If EITHER nx or ny is GREATER than the minimum values,
; 		then a new (finer) cartesian mesh is made and a bi-cubic spline
; 		routine is used to tabulate psi on the new mesh.
; 		BEWARE the bi-cubic spline routine is time consuming.
; KEYWORDS --
;     	do_plot -- if this keyword is set then the routine produces 
;          	graphical output, as specified below.
;     	do_one -- if this keyword is set then the routine calculates
;          	only one contour, corresponding to psi = psi2.
; OUTPUTS --
;     	psi = values of psi where contours have been determined,
;           	where psi is the NORMALIZED poloidal flux, so that
;           	psi = 0 in the center and psi = 1 at the boundary.
;           	THE DEFAULT is 30 surfaces, equally spaced in psi 
;           	from psi1 to psi2. (Note, 30 is the maximum number 
;           	of contours that the IDL contour routine will make 
;           	in a single call. If psi1 is too small, the routine 
;           	may not be able to find the corresponding contour,  
;           	in which case only 29 contours will be returned.)
;           	OTHERWISE, if the DO_ONE keyword is set, then
;           	only one contour (for psi = psi2) is calculated.
;     	thetag = vector defining theta grid (equally spaced 
;              	from 0 to 2*pi, not including 2*pi).
;     	rg = array giving the minor radius r(theta,psi) at each
;          	theta grid point on a given psi surface -- in METERS.
;     	b = B-field components in polar (phi,theta,r) coordinates:
;         	b(*,*,0) = B_phi, b(*,*,1) = B_theta, b(*,*,2) = B_radial,
;         	where each component, B_phi, etc. is an array of the same 
;         	dimensions as rg.  The units of b are TESLA.
;
; GRAPHICAL OUTPUT -- two plots --
;     	(1) plot of qbar versus psi, where qbar is the theta average 
;         	of local safety factor, qlocal. (Note qbar = q_MHD).
;         	The '+' symbols mark the values of qbar given 
;         	by the EFIT calculation (as stored in g.qpsi).
;     	(2) polar plot of flux surfaces, in the form r = f(theta,psi),
;         	where f(theta,psi) is stored in the array rg.
; NOTE ON COORDINATE SYSTEMS -- 
; 	If (phi,R,Z) denotes the usual cylindrical coordinates then
; 	theta and r are defined as follows: R = R_o + r*cos(theta) and 
; 	Z = Z_o + r*sin(theta). Here R_o and Z_o are the position of the
; 	magnetic axis, which is Z_o = g.zmaxis and R_o = g.rmaxis.
; 	The flux surfaces are specified by r = r(theta,psi).
;
;************************************************************************

pro contour_psi, g,npts,psi1,psi2,nx,ny, $
                 psi,thetag,rg,b, $
                 do_plot=do_plot,do_one=do_one
compile_opt defint32,strictarr,strictarrsubs
;
; Define psixy on (x,y) cartesian grid 
; psixy = the poloidal flux function on (x,y) grid,
; Here psixy(i,j) is psi at x = x(i) and y = y(j).
;

x = linspace(g.rgrid1,g.rgrid1 + g.xdim,g.mw)
y = linspace(g.zmid - g.zdim/2.,g.zmid + g.zdim/2.,g.mh)
psixy = g.psirz[0:g.mw-1,0:g.mh-1]
x = temporary(x) - g.rmaxis
y = temporary(y) - g.zmaxis

if (nx gt g.mw) or (ny gt g.mh) then begin
  ;
  ; use bi-cubic spline to increase the number of (x,y) grid 
  ; points where psixy is specified.
  ;
  x0 = x
  y0 = y 
  x = linspace(min(x0),max(x0),nx)
  y = linspace(min(y0),max(y0),ny)
  psixy = bispline_grid(psixy,x0,y0,x,y)
endif

;
; calculate the psi contours in polar (r,theta) coordinates,
; using the idl contour routine.

dpsi = g.ssibry - g.ssimag

if keyword_set(do_one) then $ 
  psi = g.ssimag + psi2*dpsi $
else $
  psi = linspace(g.ssimag + psi1*dpsi,g.ssimag + psi2*dpsi,30)

;
; These are the surfaces for which contours will be found.
; ** psi must be a vector of length less than or equal to 30 **
;

thetag = 2.*!pi*findgen(npts)/float(npts)
rg = fltarr(npts,n_elements(psi))

;
; Because the system parameters used by the convert_coord routine
; are set by the most recent graphics output, we must rescale them
; (to order unity) to avoid possible floating point errors when 
; making the coordinate transformation. Note, if the overplot
; keyword is not set when calling the contour routine, then the
; rescaling is done automatically, but the current graphics window
; is cleared, which should generally be avoided.

;  tbt: But /overplot was causing other problems - I took it out.
;

;Print, 'Contour_psi: 10/31/96  New graphic reset correction.'

old_x = !x
old_y = !y
;!x.s = [0.,1.]
;!y.s = [0.,1.]

;  Reset all scaling items to original default vaules.

   !x.TITLE           = ' '      ; STRING    ''
   !x.TYPE            = 0        ; LONG                 0
   !x.STYLE           = 0        ; LONG                 0
   !x.TICKS           = 0        ; LONG                 0
   !x.TICKLEN         = 0.       ; FLOAT           0.00000
   !x.THICK           = 0.       ; FLOAT           0.00000
   !x.RANGE           = 0.       ; FLOAT     Array(2)
   !x.CRANGE          = 0.       ; FLOAT     Array(2)
   !x.S               = [ 0.,1.] ; FLOAT     Array(2)
   !x.MARGIN          = [10.,3.] ; FLOAT     Array(2)
   !x.OMARGIN         = 0.       ; FLOAT     Array(2)
   !x.WINDOW          = 0.       ; FLOAT     Array(2)
   !x.REGION          = 0.       ; FLOAT     Array(2)
   !x.CHARSIZE        = 0.       ; FLOAT           0.00000
   !x.MINOR           = 0        ; LONG                 0
   !x.TICKV           = 0.       ; FLOAT     Array(30)
   !x.TICKNAME        = ' '      ; STRING    Array(30)
   !x.GRIDSTYLE       = 0        ; LONG                 0
   !x.TICKFORMAT      = ' '      ; STRING    ''


   !y.TITLE           = ' '      ; STRING    ''
   !y.TYPE            = 0        ; LONG                 0
   !y.STYLE           = 0        ; LONG                 0
   !y.TICKS           = 0        ; LONG                 0
   !y.TICKLEN         = 0.       ; FLOAT           0.00000
   !y.THICK           = 0.       ; FLOAT           0.00000
   !y.RANGE           = 0.       ; FLOAT     Array(2)
   !y.CRANGE          = 0.       ; FLOAT     Array(2)
   !y.S               = [ 0.,1.] ; FLOAT     Array(2)
   !y.MARGIN          = [10.,3.] ; FLOAT     Array(2)
   !y.OMARGIN         = 0.       ; FLOAT     Array(2)
   !y.WINDOW          = 0.       ; FLOAT     Array(2)
   !y.REGION          = 0.       ; FLOAT     Array(2)
   !y.CHARSIZE        = 0.       ; FLOAT           0.00000
   !y.MINOR           = 0        ; LONG                 0
   !y.TICKV           = 0.       ; FLOAT     Array(30)
   !y.TICKNAME        = ' '      ; STRING    Array(30)
   !y.GRIDSTYLE       = 0        ; LONG                 0
   !y.TICKFORMAT      = ' '      ; STRING    ''

; Jeff Schachter 1998.08.06: NOTE THAT ![X,Y].TICKFORMAT = ' ' PREVENTS
;                            THE USER FROM PLOTTING!  THIS CAUSES IDL TO 
;			     SEARCH FOR A PROCEDURE NAMED ' ' TO PROVIDE
;			     THE TICKLABELS.  KUPFER SHOULD HAVE DONE 
;			     ![X,Y].TICKFORMAT = ''

psi = psi[sort(psi)]
contour, psixy, x, y, levels=psi, path_info=header,path_xy=path_xy, /overplot

for i=0,n_elements(header)-1 do begin
   if (header[i].type eq 1) and (header[i].n ge 6) then begin
    xyarr = path_xy[*,header[i].offset:header[i].offset+header[i].n-1]


;  Jeff Schachter 1998.08.06:
;	There is no need to make a call to CONVERT_COORD here!!!
;	Since there is no plot defined, and since Kupfer set !X.S = [0,1] = !Y.S
;	above, data and normalized coordinates are identical.
;	The array xy_t[0,*] is IDENTICAL to xyarr[0,*] and likewise for xy_t[1,*]
;	Since starting in IDL 5.1 calls to CONVERT_COORD are not allowed without
;	a coordinate system defined by creating a plot (this could still be a bug
;	in IDL), and since the CONVERT_COORD is unnecessary anyway, I have commented
;	it out and replaced it with a direct equation of the important arrays.
;    xy_t = convert_coord(xyarr(0,*),xyarr(1,*),/normal,/to_data)
;    x_c = xy_t(0,*)
;    y_c = xy_t(1,*)

	x_c = xyarr[0,*]	; Jeff Schachter 1998.08.06
	y_c = xyarr[1,*]	; Jeff Schachter 1998.08.06

    ;
    ; -- Note -- to make an (x,y) plot of the original contours, 
    ; remove the overplot keyword from the contour command and 
    ; insert here the following command: 
    ; oplot, x_c, y_c
    if (min(x_c)*max(x_c) lt 0.) and (min(y_c)*max(y_c) lt 0.) then begin
      ;
      ; -- interpolate contour into (r,theta) coordinates --
      ;
      r = sqrt(x_c^2 + y_c^2)
      theta = acos(x_c/r)
      z = where(y_c lt 0.)
      theta[z] = 2.*!pi - temporary(theta[z])
      z = uniq(theta,sort(theta))
      theta = temporary(theta[z])
      r = temporary(r[z])
      n = n_elements(r)
      r_s = fltarr(n+2)
      theta_s = fltarr(n+2)
      theta_s[1:n] = theta
      r_s[1:n] = r
      r_s[0] = r[n-1]
      theta_s[0] = theta[n-1] - 2.*!pi
      r_s[n+1] = r[0]
      theta_s[n+1] = theta[0] + 2.*!pi
      rg[*,header[i].level] = spline(theta_s,r_s,thetag)
    endif
  endif
endfor

;endwhile
;
;free_lun, unit   ;close file

!x = old_x
!y = old_y

;!x.s = old_x     ;reset previous graphics scaling
;!y.s = old_y

i = where(rg[0,*] ne 0.)   ; remove zeroes

if i[0] eq -1 then begin
  ;
  ; -- ERROR condition -- 
  ;
  print, 'Contour_psi: Error - could not find any contours.'
  return
endif

;rg  =  temporary(rg(*,i))
 rg  =  temporary(rg[*,i])

; This fixes a problem of Degrassi who had a function psi.pro
;psi =  temporary(psi(i))
 psi =  temporary(psi[i])

psi = (temporary(psi) - g.ssimag)/dpsi ; normalize psi

;
; -- calculate B-field components in (phi,theta,r) basis --
;

metric, psixy,x,y,thetag,rg,hxy_g,gxy_g
numpts = g.mw
psi_x  = linspace(0.,1.,numpts)
fpsi_x = g.fpol[0:(numpts-1)]
fpsi   = spline(psi_x,fpsi_x,psi)

s1  = n_elements(rg[*,0])
s2  = n_elements(rg[0,*])
x_g = fltarr(s1,s2)
for n=0,s2-1 do x_g[*,n] = rg[*,n]*cos(thetag)

brad   = gxy_g/(x_g + g.rmaxis)/rg
btheta = sqrt(hxy_g^2 - (gxy_g/rg)^2)/(x_g + g.rmaxis)
btor   = fltarr(s1,s2)
for n=0,s2-1 do btor[*,n] = abs(fpsi[n])/(x_g[*,n] + g.rmaxis)

b = fltarr(s1,s2,3)
b[*,*,0] = btor
b[*,*,1] = btheta
b[*,*,2] = brad

if keyword_set(do_plot) then begin 
  ;
  ; -- BEGIN GRAPHICAL OUTPUT --
  ;
  ; -- calculate local safety factor --
  ;
  qlocal = btor*rg/btheta/(x_g + g.rmaxis)
  ;
  ; -- compute qbar (the theta average of qlocal) -- 
  ;
  qbar = fltarr(s2)
  for n=0,s2-1 do qbar[n] = total(qlocal[*,n])/float(s1)
  ;
  ; -- plot qbar versus psi, comparing to the pre-computed results
  ;     stored in g.qpsi (i.e. output from EFIT) -- 
  ;
  !p.multi = [0,2,1] ; makes two plots across page.
  plot, [0.,1.], [0.,max(qbar)], /nodata, xtitle='psi', ytitle='q'
  oplot, psi, qbar, psym=1 
  oplot, psi_x[0:[numpts-2]], g.qpsi[0:[numpts-2]]
  ;
  ; -- plot the flux surfaces in polar coordinates --
  ;
  n = s2-1
  plot, /polar,[rg[*,n],rg[0,n]],[thetag,2.*!pi]
  for n=0,s2-2 do oplot,/polar,[rg[*,n],rg(0,n)],[thetag,2.*!pi]
  !p.multi = 0 ; return to single-plot-per-page mode.
  ;
  ; -- END GRAPHICAL OUTPUT --
  ;
endif

return
end
