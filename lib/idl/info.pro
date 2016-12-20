PRO info, str
    ;+#info
    ;+Print a informational message
    ;+***
    ;+##Arguments
    ;+     **str**: message
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> info, "This is an informative message"
    ;+```
    print, colored('INFO: '+str,c='b')
END
