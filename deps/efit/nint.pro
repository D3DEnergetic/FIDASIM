FUNCTION NINT,r_val
;
; Convert to nearest integer value.  Handles negative numbers
;
; 06-Nov-00 Sped up [BD]
; 02-Nov-00 Fix for array of 1 [BD]
; 07-Sep-89 BD Written
;
info = SIZE( r_val )
if info[0] EQ 0 THEN BEGIN
   IF r_val GE 0 THEN RETURN,LONG(r_val+0.5)
   RETURN,LONG(r_val-0.5)
ENDIF ELSE BEGIN
   nints =  LONG(r_val+0.5)
   
   negInds = WHERE( r_val LT 0.0, nNeg )

   IF nNeg GT 0 THEN nints( negInds ) = LONG(r_val( negInds )-0.5)

RETURN,nints
ENDELSE
END
