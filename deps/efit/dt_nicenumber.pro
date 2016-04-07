;+-------------------------------------------------------------------
;  NAME:
; 	dt_nicenumber
;
;  PURPOSE:
; 	Makes nice delta time numbers, by rounding to 3 significant digits
;	Can be used if you want to know if you have a constant timebase.
;
; CATEGORY:
;       Math, Graphics
;
; CALLING SEQUENCE:
;        niceDts = dt_nicenumber( dts )
;
; INPUTS:
;	dts - an array of numbers
;
; KEYWORD PARAMETERS:
;    Optional Inputs:
;	nSignificantDigits - number of significant digits to round to
;				(defaults to 3)
;    Optional Outputs:
;       nDecimals - number of significant decimal places in returned value
;  EXAMPLE:
;	IDL> dt=[2.0003,2.15,.2003,.20009,22.0004,.250002,.5001,100000.1]
;	IDL> print, dt_nicenumber(dt) 
;	   2.00000    2.15000   0.200000   0.200000    22.0000   0.250000
;	   0.500000    100000.
;
;	IDL> print,  dt_nicenumber( -1.876e6 ) 
; NOTES:
;	see nicenumber.pro to round to numbers like 1, 2, 2.5, 5, etc.,
;	which might be time intervals for older digitizer rates.
;  HISTORY:
;	13-May-2011 fixed bug
;	20-Apr-2011 handle arrays
;	05-Nov-2008 corrected for numbers like -1.876e6
;	26-Sep-2005 fixed nToShift for nums less than zero (like 0.0016)
;	26-Jul-04 work OK for 0
;  	03-Nov-00 Rewrote [BD]
;--------------------------------------------------------------------
function dt_nicenumber, dt_in, nSignificantDigits=nSignificantDigits, $
                        ndecimals=nDecimals	;nToShift

on_error, 2

; for delta t rounding
;  (can assume to be positive)

if n_elements(nSignificantDigits) EQ 0 THEN nSignificantDigits = 3
nToShift = nSignificantDigits	; default if dt_in eq 0.0

if dt_in[0] eq 0.0 and n_elements(dt_in) eq 1 then begin
   nDecimals = 3
   return, dt_in	; in case used like this
endif

nice =  dt_in

for i=0, n_elements( dt_in )-1 do begin

   if nice[i] ne 0.0 then begin

      posDT = abs( dt_in[i] )
      if dt_in[i] lt 0 then sign = -1 else sign = 1

      exponent =  ALOG10( posDT )
      dt = DOUBLE( posDT )
      varInfo = SIZE( posDT )

      ;;;negInds = WHERE( exponent LE 3, nNeg )
      ;;;posInds = WHERE( exponent GT 4, nPos )
      negInds = WHERE( exponent LT 0, nNeg )
      posInds = WHERE( exponent GE 0, nPos )

      nice[i] = DOUBLE( NINT(dt) )	; to handle the exponent of 0

      if nNeg GT 0 THEN BEGIN
	 nToShift = FIX( nSignificantDigits - exponent )
	 nice[i] = ( NINT(dt * 10D^nToShift) ) / 10D^nToShift

	 if arg_present( nDecimals ) then begin	; probably slow, so only if needed
	    str = string( nice[i])
	    for ipart=0,n_elements( str )-1 do begin
	       temp=str[ipart]
	       for ic=strlen(temp),1,-1 do begin
        	  if  (strmid(temp, ic-1, 1 ) ne '0') then goto, done
	       endfor
      done:
	       trimmed = strmid( temp, 0, ic)
	       nToShift = strlen( strtrim(trimmed,2) ) -2
	    endfor
	 endif

      ENDIF

      if nPos GT 0 THEN BEGIN
	 nToShift = FIX( exponent) - nSignificantDigits +1
	 nice[i] = DOUBLE( NINT(dt/ 10.^nToShift) ) * 10.^nToShift
      ENDIF

      if arg_present( nDecimals ) then begin
	 nDecimals = nToShift
	 if n_elements( nDecimals ) eq 1 then nDecimals = nDecimals[0]
      endif

      if sign eq -1 then nice[i] = nice[i]*sign
   end


end

if varInfo[1] EQ 5 THEN return, nice ELSE return, FLOAT(nice)

end
