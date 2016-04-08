PRO error, str
    ;+#error
    ;+Print a error message
    ;+***
    ;+##Arguments
    ;+     **str**: message
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> error, "=("
    ;+```
    print, colored('ERROR: '+str,c='r')
END
