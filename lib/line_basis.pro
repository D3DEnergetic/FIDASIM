FUNCTION line_basis, r0, v0, inv_basis=inv_basis
    ;+#line_basis
    ;+Calculates basis from a line with +x in the direction of line
    ;+***
    ;+##Arguments
    ;+    **r0**: Starting point of line [cm]
    ;+
    ;+    **v0**: Direction of line
    ;+
    ;+##Keyword Arguments
    ;+    **inv_basis**: Set this to a named variable that recieves the inverse basis
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> basis = line_basis([0,0,0],[0,-1,0])
    ;+IDL> x = basis##[1,1,0] ;Transforms a point in line-space ([1,1,0]) to real space
    ;+IDL> print,x
    ;+    [1, -1, 0]
    ;+```

    rf = r0 + v0
    dis = sqrt(total(v0^2.0))
    beta = asin((r0[2] - rf[2])/dis)
    alpha = atan((rf[1] - r0[1]),(rf[0]-rf[0]))

    R = tb_zyx(alpha,beta,0.0)
    inv_basis = transpose(R)
    return, R

END
