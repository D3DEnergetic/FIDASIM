;***********************************************************************
;
;  4dlib/EFITLIB/kupfer/idl_geqdsk/fluxmap.pro
;
;  created:
;
;  modified:
;    4/02/2002  TBT - put incheck for ni: i = where(y_bdry lt 0., ni)
;    6/23/97    Ted Terpstra - put in check for numpts = 0
;    8-9-96     K. Greene -- Add use of function TEMPORARY
;
;
;**********************************************************************
;
; -- Calculates the radial flux coordinate (rho) and tabulates 
;    the poloidal flux function (psi) on an equally spaced rho grid.
;
; -- INPUT --
; -- g = structure obtained from readg containing the GEQDSK.
; -- OUTPUT --
; -- output = data structure containing the following fields:
; -- shot = shot number (g.shot)
; -- time = time in msec (g.time)
; -- rho = values of rho IN METERS at 101 equally spaced 
;          grid points, from the magnetic axis to the boundary.
;          Note rho = sqrt[(Toroidal Flux)/(pi*bcentr)], where
;          bcentr is the toroidal field on axis -- from the 
;          GEQDSK -- we take bcentr = abs(g.bcentr).
; -- psi = corresponding values of poloidal flux function (in MKS units).
;       psi = g.ssimag at the magnetic axis,
;       psi = g.ssibry at the boundary (last closed flux surface).
; -- bcentr = abs(g.bcentr) in Tesla, as used to define rho.
;
; The map from psi to rho is determined by the relation
;      d(Toroidal Flux)/d(psi) = 2*pi*q(psi)
; where q(psi) is the safety factor. Values of q(psi)
; are taken from the GEQDSK and the toroidal flux function
; is computed by direct integration. However, the boundary
; point (in psi) is excluded because the local safety 
; factor diverges on the separatrix (at the x-point) and
; the corresponding value of q(psi) in the GEQDKS is not
; reliable. Thus we compute the toroidal flux inside the 
; boundary by direct integration of B_phi over the the 
; area enclosed by the boundary surface, whose contour 
; is obtained directly from the GEQDSK.
;
;***********************************************************************

function fluxmap, g

numpts = g.mw

If (numpts Le 1) Then Begin
  Print, 'Fluxmap: Errot: g.mw = numpts is NOT > 1. Return ', numpts
  Return, {ierr:1}
EndIf

dpsi = g.ssibry - g.ssimag
psi_eqdsk = findgen(numpts-1)/float(numpts-1)
q_eqdsk = g.qpsi(0:(numpts-2))   

;
; -- q_eqdsk = q(psi) evaluated at the points psi = psi_eqdsk,
;              as given in the GEQDSK.
; It is assumed that the GEQDSK evaluates flux functions
; (such as q) on a grid of NUMPTS equally spaced points in psi from
; g.ssimag to g.ssibry. Here we have excluded the last point 
; (at the boundary) in our forming of psi_eqdsk and q_eqdsk.
; Also we are working in NORMALIZED poloidal flux, where
; psi=0 (at g.ssimag) and psi=1 (at g.ssibdry).
;
; -- interpolate q_eqdsk to a finer psi grid --
;

nint = 100 ; number of points on fine psi grid.
psi = psi_eqdsk(numpts-2)*findgen(nint)/float(nint-1)
q = spline(psi_eqdsk,q_eqdsk,psi) 

;
; -- integrate q(psi) to get the Toroidal Flux --
; Here we use Simpson's rule.
;

flux = 2.*qsimpv(psi,q)*abs(dpsi)   

;
; NOTE flux = (Toroidal Flux)/pi.
; Also, the toroidal flux is not evaluated at the boundary
; surface, but only to some nearby interior surface, 
; where psi = psi_eqdsk(numpts-2).
;
; -- calculate the toroidal flux for boundary surface --
; This is done in 4 STEPS, as given below.
;
; STEP 1 -- calculate the contour of an interior surface --
;

psi_c = psi_eqdsk(numpts-2) 	; interior surface
pts = 101  			; number of theta points (** an odd number **)

contour_psi, g,pts,psi_c,psi_c,g.mw,g.mh,psi_c,thetag,r_c,/do_one

;
; STEP 2 -- define boundary surface from GEQDSK --
;
if g.nbdry eq 0.0 then g.nbdry = n_elements( g.bdry(0,*) )
x_bdry = g.bdry(0,0:(g.nbdry-1)) - g.rmaxis
y_bdry = g.bdry(1,0:(g.nbdry-1)) - g.zmaxis

;
; STEP 3 -- interpolate boundary into (r,theta) coordinates --
;

r = sqrt(x_bdry^2 + y_bdry^2)
theta = acos(x_bdry/r)

; Check put in for ni - 20020402 tbt
i = where(y_bdry lt 0., ni)

If (ni Gt 0) Then theta(i) = 2.*!pi - temporary(theta(i))

i = uniq(theta,sort(theta))
theta = temporary(theta(i))
r = temporary(r(i))
n = n_elements(r)
r_s = fltarr(n+2)
theta_s = fltarr(n+2)
theta_s(1:n) = theta
r_s(1:n) = r
r_s(0) = r(n-1)
theta_s(0) = theta(n-1) - 2.*!pi
r_s(n+1) = r(0)
theta_s(n+1) = theta(0) + 2.*!pi
r_bdry = spline(theta_s,r_s,thetag)

;
; STEP 4 -- calculate the toroidal flux between the interior
;           surface and the boundary surface --
; Here the toroidal flux is calculated by direct 
; integration in polar (r,theta) coordinates.
; It is assumed that the boundary surface and the 
; interior surface are close, in the sense that 
; f(psi) = R*B_phi is nearly constant. Once f(psi)
; is taken outside of the integrand, the radial 
; integration can be done analytically. The theta 
; integration is then done simply as a summation. 
;

eps = r_c/g.rmaxis
y1 = (eps - alog(1.+eps*cos(thetag))/cos(thetag))/cos(thetag)
eps = r_bdry/g.rmaxis
y2 = (eps - alog(1.+eps*cos(thetag))/cos(thetag))/cos(thetag)
fpsi = 0.5*(g.fpol(numpts-2) + g.fpol(numpts-1))
r_integral = abs(fpsi)*g.rmaxis*(y2 - y1)
dflux = 2.*total(r_integral)/float(pts)
flux_bdry = flux(nint-1) + dflux

;
; NOTE flux_bdry = (Toroidal Flux inside boundary)/pi.
;
; -- define the equally spaced rho grid (in meters) -- 
;

bcentr = abs(g.bcentr)
rho = sqrt(flux_bdry/bcentr)*findgen(101)/100.

;
; -- tabulate psi and q on the rho grid --
;

flux_new = bcentr*rho^2
psi = spline([flux,flux_bdry],[psi,1.],flux_new)

;
; -- convert psi to un-normalized form --
;

psi = g.ssimag + dpsi*temporary(psi)

output = {fluxmap, shot:g.shot, time:g.time, $
          psi:psi, rho:rho, bcentr:bcentr}

return, output
end
