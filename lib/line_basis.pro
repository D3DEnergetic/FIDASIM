FUNCTION line_basis, r0, v0, inv_basis=inv_basis

    rf = r0 + v0
    dis = sqrt(total(v0^2.0))
    beta = asin((r0[2] - rf[2])/dis)
    alpha = atan((rf[1] - r0[1]),(rf[0]-rf[0]))

    R = tb_zyx(alpha,beta,0.0)
    inv_basis = transpose(R)
    return, R

END
