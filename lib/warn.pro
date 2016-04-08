PRO warn, str
    ;+##`warn, str`
    ;+Print a warning message
    ;+###Arguments
    ;+     **str**: message
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> warn, "This may be a problem"
    ;+```
    print, colored('WARNING: '+str,c='y')
END
