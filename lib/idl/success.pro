PRO success, str
    ;+##`success, str`
    ;+Print a success message
    ;+###Arguments
    ;+     **str**: message
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> success, "Yay!!!"
    ;+```
    print, colored('SUCCESS: '+str,c='g')
END
