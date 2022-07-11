FUNCTION cocos_struct,index
    ;+#cocos_struct
    ;+Creates a COCOS structure from a given index
    ;+***
    ;+##Arguments
    ;+    **index**: COCOS index
    ;+
    ;+##Return Value
    ;+COCOS structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> cocos = cocos_struct(5)
    ;+```
    
    exp_Bp = index lt 10 ? 0 : 1
    sigma_Bp = where(index eq [1,2,5,6,11,12,15,16]) ne -1 ? 1 : -1
    sigma_RphZ = index mod 2 ne 0 ? 1 : -1
    sigma_rhothph = where(index eq [1,2,7,8,11,12,17,18]) ne -1 ? 1 : -1
    sign_q = sigma_rhothph
    sign_pprime = -1 * sigma_Bp

    cocos = {cocos:index, exp_Bp:exp_Bp, sigma_Bp:sigma_bp, sigma_RphZ:sigma_RphZ, sigma_rhothph:sigma_rhothph, sign_q:sign_q, sign_pprime:sign_pprime}
    return,cocos
END
