;************************************************************************
;
;  4dlib/EFITLIB/kupfer/idl_geqdsk/metric.pro
;
;  created:
;
;  modified:
;
;
;***********************************************************************
;
;  Calculate Metric Coefficients in (theta,psi) coordinates.
;
; -- INPUTS -- 
;    psixy = the poloidal flux function on (x,y) grid (in MKS units).
;    x and y = vectors defining (x,y) grid -- in METERS.
;    thetag = vector defining theta grid.
;    rg = array giving the minor radius r(theta,psi) at each
;         theta grid point on a given psi surface -- in METERS.
; -- OUTPUTS --
;     hxy_g = hxy interpolated onto the (theta,psi) grid,
;        where hxy = |grad(psi)|.
;     gxy_g = gxy interpolated onto the (theta,psi) grid,
;        where gxy = r^2[grad(psi)*grad(theta)].
;
;**********************************************************************

pro metric, psixy,x,y,thetag,rg, hxy_g, gxy_g

;
; -- calculate hxy and gxy on (x,y) grid --
;

s = size(psixy)
psixy_x = fltarr(s(1),s(2))
psixy_y = fltarr(s(1),s(2))
gxy = fltarr(s(1),s(2))
for n=0,s(1)-1 do psixy_y(n,*) = deriv(y,psixy(n,*))
for n=0,s(2)-1 do psixy_x(*,n) = deriv(x,psixy(*,n))
for n=0,s(1)-1 do gxy(n,*) = x(n)*psixy_y(n,*) - y*psixy_x(n,*)
hxy = sqrt(psixy_x^2 + psixy_y^2) 

;
; -- interpolate hxy and gxy onto (theta,psi) grid --
;

s1 = n_elements(rg(*,0))
s2 = n_elements(rg(0,*))
gxy_g = fltarr(s1,s2)
hxy_g = fltarr(s1,s2)
x_g = fltarr(s1,s2)
y_g = fltarr(s1,s2)
for n=0,s2-1 do x_g(*,n) = rg(*,n)*cos(thetag)
for n=0,s2-1 do y_g(*,n) = rg(*,n)*sin(thetag)
dx = x(1) - x(0)
dy = y(1) - y(0)
gxy_g = interpolate(gxy,(x_g-x(0))/dx,(y_g-y(0))/dy)
hxy_g = interpolate(hxy,(x_g-x(0))/dx,(y_g-y(0))/dy)

return
end
