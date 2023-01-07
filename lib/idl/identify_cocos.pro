FUNCTION identify_cocos,g,ccw_phi=ccw_phi,exp_Bp=exp_Bp
    ;According to IDl documentaion, setting /BOOLEAN_KEYWORD=0 in the function definition set the default value to 1?
    ;+#identify_cocos
    ;+Identifies the cocos index of a GEQDSK dictionary object
    ;+Reference:
    ;+    O. Sauter and S. Ty. Medvedev, Tokamak Coordinate Conventions: COCOS, 
    ;+    Computer Physics Communications 184, 293 (2013).
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK structure
    ;+
    ;+##Keyword Arguments
    ;+    **ccw_phi**: Toroidal direction from top view, True if CCW, False if CW
    ;+
    ;+    **exp_bp**: 0 if poloidal flux divided by 2 pi, 1 if using effective poloidal flux
    ;+
    ;+##Return Value
    ;+COCOS structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> g = readg(filename)
    ;+IDL> cocos = identify_cocos(g)
    ;+```
   
    ;Sauter, eq. 22
    sigma_Bp_in = -1 * signum(g.pprime[0] * g.current)
    sigma_RphZ_in = ~keyword_set(ccw_phi) ? 1 : -1
    sigma_rhothph_in = signum(g.qpsi[0] * g.current * g.bcentr)

    sigmas = [sigma_Bp_in, sigma_RphZ_in, sigma_rhothph_in]

    if array_equal(sigmas, [1,1,1]) then begin
        index = 1
    endif else if array_equal(sigmas, [1,-1,1]) then begin
        index = 2
    endif else if array_equal(sigmas, [-1,1,-1]) then begin
        index = 3
    endif else if array_equal(sigmas, [-1,-1,-1]) then begin
        index = 4
    endif else if array_equal(sigmas, [1,1,-1]) then begin
        index = 5
    endif else if array_equal(sigmas, [1,-1,-1]) then begin
        index = 6
    endif else if array_equal(sigmas, [-1,1,1]) then begin
        index = 7
    endif else if array_equal(sigmas, [-1,-1,1]) then begin
        index = 8
    endif else begin
        index = FIDASIM_default_COCOS
    endelse

    cocos = ~keyword_set(exp_Bp) ? cocos_struct(index) : cocos_struct(index + 10)
    return,cocos
END
