;***********************************************************************
;
;  4dlib/EFITLIB/kupfer/idl_math/linspace.pro
;
;  created:
;
;  modified:
;
;
;***********************************************************************

function linspace, a,b,n

compile_opt defint32,strictarr,strictarrsubs

;
; create a vector x, with n evenly spaced steps from a to b
;

if n ge 2 then $
  x = a + (b-a)*findgen(n)/float(n-1) $
else x = (a + b)/2.

return,x
end
