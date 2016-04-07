;-------------------------------------------------------------
;
; NAME:
;       DELCHR
; PURPOSE:
;       Delete all occurrences of a character from a text string.
; CATEGORY:
;	Programming
; CALLING SEQUENCE:
;       new = delchr(old, char)
; INPUTS:
;       old = original text string.     in
;       char = character to delete.     in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       new = resulting string.         out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner.  5 Jul, 1988.
;       Johns Hopkins Applied Physics Lab.
;       RES 11 Sep, 1989 --- converted to SUN.
;       R. Sterner, 27 Jan, 1993 --- dropped reference to array.
;
; Copyright (C) 1988, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
        FUNCTION DELCHR, OLD, C, help=hlp
 
        if (n_params(0) lt 2) or keyword_set(hlp) then begin
          print,' Delete all occurrences of a character from a text string.'
          print,' new = delchr(old, char)
          print,'   old = original text string.     in'
          print,'   char = character to delete.     in'
          print,'   new = resulting string.         out'
          return, -1
        endif
 
        B = BYTE(OLD)                      ; convert string to a byte array.
        CB = BYTE(C)                       ; convert char to byte.
        w = where(b ne cb(0))
        if w(0) eq -1 then return, ''      ; Nothing left.
        return, string(b(w))               ; Return new string.
        END
