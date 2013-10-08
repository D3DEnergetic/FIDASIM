;**************************************************************************
;
;  4dlib/EFITLIB/kupfer/idl_geqdsk/rho_rz.pro
;
;  created:
;
;  modified:
;    4/02/2002  tbt - put in check for i >0: theta(i) = 2.*!pi - theta(i)
;    7/27/00    tbt - bombing with bad index in gaprofiles/density for
;               72242/3400. Caused by readg structure having 1
;               dimension for g.psirz for 72242.
;             
;    9/24/96    Ted Terpstra - idl was sometimes going into
;                 an infinite run state when r(i)=0 for some i
;                 in acos(x_pts/r). Put in check to stop /0.0
;    8-09-96     K. Greene -- Add use of function TEMPORARY
;    2/05/96    TBT - Put in check for # of r_pts = # z_pts
;    6-28-95	K. Greene -- changed names of numerical recipe
;		references (IDL 4.0)
;			nr_spline   -->   spl_init
;			nr_splint   -->   spl_interp
;
;
;************************************************************************
;
; -- This function calculates rho at specified points in
;    the (R,Z) plane, where R is the major radius and Z is
;    the height above the midplane.
;
; -- INPUT --
; -- g = structure obtained from readg containing the GEQDSK.
; -- r_pts = the R locations for a set of points (specified below)
; -- z_pts = the Z locations of a set of points (specified below)
;      r_pts(i) = R location of the i-th point (in METERS).
;      z_pts(i) = Z location of the i-th point (in METERS).
;
; -- KEYWORDS --
; -- do_linear -- if this keyword is present and non-zero, then
;        the routine uses bi-linear interpolation, instead of the
;        more time consuming bi-cubic spline.
;
; -- norm -- if set and non zero, the normalized rho is returned.
;        In addition, values outside the separatrix are in terms
; 	 of the square root of normalized flux.  (CBF)
;
; -- OUTPUT --
; -- rho_pts = the rho value in METERS corresponding to the 
;      specified (R,Z) points (as given above) such that
;      rho_pts(i) = rho at the i-th point.
; -- psi_pts = psi values corresponding to rho_pts (in MKS units).
; -- index = indices of points, [r_pts(index),z_pts(index)], 
;      which are outside the last closed flux surface. By convention,
;      we set rho_pts(index) = rho_at_boundary for all such points.
;      If no points are outside boundary, then index = -1.      
; -- rhobnd = the rho value at the boundary. (CBF)
;
; -- define coordinates (as above) --
;
;***********************************************************************

function rho_rz, g,r_pts, z_pts, psi_pts, index, $
                 do_linear=do_linear,            $
		 norm = norm, rhobnd = rhobnd

compile_opt defint32,strictarr,strictarrsubs

if n_elements(norm) le 0 then norm = 1

;Print, ' '
;Print, 'Rho_rz: norm = ', norm

nx_pts = N_Elements(r_pts)
ny_pts = N_Elements(z_pts)
If (nx_pts Ne ny_pts) Then Begin
   Print, 'RHO_RZ:  ERROR - # x_pts NE # y_pts ', nx_pts, ny_pts
EndIf

x_pts = r_pts
y_pts = z_pts

; -- setup psi on the cartesian grid given by GEQDSK --
x = linspace(g.rgrid1,g.rgrid1 + g.xdim,g.mw)
y = linspace(g.zmid - g.zdim/2.,g.zmid + g.zdim/2.,g.mh)
psixy = g.psirz[0:g.mw-1,0:g.mh-1]

if keyword_set(do_linear) then begin 
  ;
  ; -- use bi-linear interpolation to get psi at (x,y) --
  ;
  dx = x[1] - x[0]
  dy = y[1] - y[0]
  psi_pts = interpolate(psixy,(x_pts-x[0])/dx,(y_pts-y[0])/dy)
endif else begin
  ;
  ; -- use bi-cubic spline interpolation to get psi at (x,y) --
  ;
  psi_pts = bispline(psixy,x,y,x_pts,y_pts)
endelse

;
; -- interpolation to get rho given psi --
;

map = fluxmap(g)
psi = map.psi
rho = map.rho

;rho_pts = spl_interp(psi,rho, spl_init(psi,rho),psi_pts)
rho_pts   = nr_splint(psi,rho,nr_spline(psi,rho),psi_pts)

;
; -- now find all points that are outside the last closed flux surface 
; -- and set rho_pts = max(rho) for all such points.  
;
; -- obtain boundary points from GEQDSK --
;

x_bdry = g.bdry[0,0:(g.nbdry-1)]
y_bdry = g.bdry[1,0:(g.nbdry-1)]

;
; -- express boundary points in (r,theta) coordinates --
; Here (r,theta) denotes ordinary polar coordinates,
; where the origin (r = 0) is at the magnetic axis.
; Boundary points are (r,theta) = (r_bdry,theta_bdry).
;

x_bdry = temporary(x_bdry) - g.rmaxis
y_bdry = temporary(y_bdry) - g.zmaxis
r      = sqrt(x_bdry^2 + y_bdry^2)


; tbt fix ////////////////////////////////////////

 theta = r * 0.  ; Give theta right size = 0.

 i = where(r Ne 0.0, inum)

; Print, 'Rho_rz: i for where (r.Ne. 0.0) = ', inum
 If (inum Ne 0) Then theta[i] = acos(x_bdry[i]/r[i])

;theta = acos(x_bdry/r)

; tbt fix done ///////////////////////////////////


i = where(y_bdry lt 0., in_y_bdry)

;Print, 'rho_rz: i for where y_bdry .lt. 0 = ', in_y_bdry
;;theta(i) = 2.*!pi - theta(i)  20020402 tbt - new

If (in_y_bdry Gt 0) Then theta[i] = 2.*!pi - theta[i]

;
; -- sort and make ends periodic --
;

i               = uniq(theta,sort(theta))

;Print, 'Rho_rz: N_elements(uniq(i)) = ', N_elements(i)

theta           = theta[i]
r               = r[i]
n               = n_elements(r)
r_bdry          = fltarr(n+2)
theta_bdry      = fltarr(n+2)
theta_bdry[1:n] = theta
r_bdry[1:n]     = r
r_bdry[0]       = r[n-1]
theta_bdry[0]   = theta[n-1] - 2.*!pi
r_bdry[n+1]     = r[0]
theta_bdry[n+1] = theta[0] + 2.*!pi

;
; -- express (x_pts,y_pts) in (r,theta) coordinates --
;

x_pts = temporary(x_pts) - g.rmaxis
y_pts = temporary(y_pts) - g.zmaxis
r     = sqrt(x_pts^2 + y_pts^2)


; tbt Fix //////////////////////////////////////////////////////

theta_pts = r * 0.   ; give theta_pts the size of r and = 0.

i = where(r Ne 0.0, inum)

;Print, 'Rho_rz: inum for where r Ne 0.0 = ', inum

If (inum Ne 0) Then theta_pts[i] = acos(x_pts[i]/r[i])

;theta_pts = acos(x_pts/r)

; tbt Fix  End  //////////////////////////////////////////////////////


i = where(y_pts lt 0., inum)

;Print, 'Rho_rz: inum for where y_pts Lt 0. = ', inum

if (inum Ne 0) then theta_pts[i] = 2.*!pi - theta_pts[i]

;
; -- interpolate to get the radial position of the boundary 
;    evaluated at theta = theta_pts --
;

r_boundary =  nr_splint( theta_bdry, r_bdry,   $
              nr_spline( theta_bdry, r_bdry),  $
              theta_pts)

; spl_interp(theta_bdry,r_bdry,spl_init(theta_bdry,r_bdry),theta_pts)
; -- points lie outside boundary if (r > r_boundary) --


index = where(r gt r_boundary, inum)

;Print, 'Rho_rz: inum for where r Gt r_boundary = ', inum

rhobnd = max(rho)

if (inum Ne 0) then rho_pts[index] = rhobnd

if norm ge 1 then begin
	rho_pts = temporary(rho_pts) / rhobnd
	if (inum Ne 0) then	rho_pts[index] = $
	   ((-psi_pts[index]+g.ssimag)/(g.ssimag-g.ssibry))^0.5
endif

;
; -- output --
;

return, rho_pts
end
