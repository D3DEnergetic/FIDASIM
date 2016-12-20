PRO error, str,halt=halt
    ;+#error
    ;+Print a error message
    ;+***
    ;+##Arguments
    ;+     **str**: message
    ;+
    ;+##Keyword Arguments
    ;+     **halt**: Halt program execution
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> error, "=("
    ;+```
    if keyword_set(halt) then begin
        message, colored(str,c='r'),level=-1
    endif else begin
        print, colored('ERROR: '+str,c='r')
    endelse
END
