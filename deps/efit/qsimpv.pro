;********************************************************************
;
;  4dlib/EFITLIB/kupfer/idl_math/qsimpv.pro
;
;  created:
;
;  modified:
;
;
;*******************************************************************
;
; compute integral of a discrete function via 
; lowest order Simpson's rule --
;
; -- input function is y(x) --
; -- output is f(x) --
;
;********************************************************************

function qsimpv, x,y

n = n_elements(x)
f = fltarr(n)
i = indgen(n-1)
df = (y(i+1)+y(i))*(x(i+1)-x(i))/2;
for j=1,(n-1) do f(j) = f(j-1) + df(j-1)

return,f
end
