;+------------------------------------------------------------
;
; NAME:
;       ISNUMBER
; PURPOSE:
;       Determine if a text string is a valid number.
; CATEGORY:
; CALLING SEQUENCE:
;       i = isnumber(txt, [x])
; INPUTS:
;       txt = text string to test.                      in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       x = optionaly returned numeric value if valid.  out
;       i = test flag:                                  out
;           0: not a number.
;           1: txt is a long integer.
;           2: txt is a float.
;           -1: first word of txt is a long integer.
;           -2: first word of txt is a float.
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;	28-Feb-2011 don't call a character when has a minus sign [BD]
;	15-Feb-2008 don't consider date string, like 02/15/2008, a number [BD]
;	29-Oct-2007 return 0 if input not defined [BD]
;	01-Dec-2006 handle structures (treat as not a number) [Bill Davis]
;       Richard Garrett, 14 June, 1992 --- fixed bug in returned float value.
;       R. Sterner, 12 Mar, 1990 --- upgraded.
;       R. Sterner.  15 Oct, 1986.
;       Johns Hopkins Applied Physics Lab.
;
; Copyright (C) 1986, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-------------------------------------------------------------
 
        FUNCTION ISNUMBER, TXT0, X, help=hlp
        
	ON_ERROR,2   
        
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
          print,' Determine if a text string is a valid number.'
          print,' i = isnumber(txt, [x])
          print,'   txt = text string to test.                      in'
          print,'   x = optionaly returned numeric value if valid.  out'
          print,'   i = test flag:                                  out'
          print,'       0: not a number.'
          print,'       1: txt is a long integer.'
          print,'       2: txt is a float.'
          print,'       -1: first word of txt is a long integer.'
          print,'       -2: first word of txt is a float.'
          return, -1
        endif
	
	if n_elements(txt0) eq 0 then return, 0
	
	if size( txt0, /type ) eq 8 then return, 0
 
        TXT = STRTRIM(TXT0[0],2)   ; trim blanks.
        X = 0                   ; define X.
 
        IF TXT EQ '' THEN RETURN, 0     ; null string not a number.
 
        SN = 1
        IF NWRDS(TXT) GT 1 THEN BEGIN   ; get first word if more than one.
          SN = -1
          TXT = GETWRD(TXT,0)
        ENDIF
          
        f_flag = 0              ; Floating flag.
        b = byte(txt)
        w = where(b eq 43, cnt)
        if cnt gt 1 then return, 0
        t = delchr(txt,'+')
	   ; why call a character when has a minus sign?
;;;        w = where(b eq 45, cnt)
;;;        if cnt gt 1 then return, 0
        t = delchr(t,'-')
        w = where(b eq 46, cnt)                 ; '.'
        if cnt gt 1 then return, 0              ; May only be 1.
        if cnt eq 1 then f_flag = 1             ; If one then floating.
        t = delchr(t,'.')
        w = where(b eq 101, cnt)                ; 'e'
        if cnt gt 1 then return, 0
        if cnt eq 1 then f_flag = 1
        t = delchr(t,'e')
        w = where(b eq 69, cnt)                 ; 'E'
        if cnt gt 1 then return, 0
        if cnt eq 1 then f_flag = 1
        t = delchr(t,'E')
        w = where(b eq 100, cnt)                ; 'd'
        if cnt gt 1 then return, 0
        if cnt eq 1 then f_flag = 1
        t = delchr(t,'d')
        w = where(b eq 68, cnt)                 ; 'D'
        if cnt gt 1 then return, 0
        if cnt eq 1 then f_flag = 1
        t = delchr(t,'D')
        if total((b eq 101)+(b eq 69)+(b eq 100)+(b eq 68)) gt 1 then return,0
        b = byte(t)
        if total((b ge 65) and (b le 122)) ne 0 then return, 0
 
        c = strmid(t,0,1)
        if (c lt '0') or (c gt '9') then return, 0  ; First char not a digit.
	
	if strpos(txt,'/') ge 0 then return, 0  ; probably a date, like 02/15/2008
 
        x = txt + 0.0                               ; Convert to a float.
        if f_flag eq 1 then return, 2*sn            ; Was floating.
        if x eq long(x) then begin
          x = long(x)
          return, sn
        endif else begin
          return, 2*sn
        endelse
 
        END
