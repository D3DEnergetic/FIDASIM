PRO convert_cocos,g,cc_in,cc_out,sigma_Ip=sigma_Ip,sigma_B0=sigma_B0,l_d=l_d,l_B=l_B,exp_mu0=exp_mu0
    ;+#convert_cocos
    ;+Converts a GEQDSK structure according to cc_in --> cc_out
    ;+Reference:
    ;+    O. Sauter and S. Yu. Medvedev, Tokamak Coordinate Conventions: COCOS,
    ;+    Computer Physics Communications 184, 293 (2013).
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK structure
    ;+
    ;+    **cc_in**: COCOS structure, input
    ;+
    ;+    **cc_out**: COCOS structure, output
    ;+
    ;+## keyword Arguments
    ;+    **sigma_Ip**: Typle of current sign, (in, out)
    ;+
    ;+    **sigma_B0**: Tuple of toroidal field sign, (in, out)
    ;+
    ;+    **l_d**: Tuple of length scale, (in, out)
    ;+
    ;+    **l_B**: Tuple of field magnitude scale, (in, out)
    ;+
    ;+    **exp_mu0**: Tuple of exponents for mu0, (in, out)
    ;+
    ;+##Return Value
    ;+GEQDSK structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> g = efit.readg(filename)
    ;+IDL> cc_in = cocos_struct(3)
    ;+IDL> cc_out = cocos_struct(1)
    ;+IDL> converted_g = convert_COCOS(g, cc_in, cc_out)
    ;+```
    print,'CONVERT_COCOS: cocos_in (',cc_in.cocos,') != cocos_out (',cc_out.cocos,') applying COCOS conversion.'
    mu0 = 4*!dpi*1e-7

    l_d_eff = l_d[1]/l_d[0]
    l_B_eff = l_B[1]/l_B[0]
    exp_mu0_eff = exp_mu0[1] - exp_mu0[0]

    exp_Bp_eff = cc_out.exp_Bp - cc_in.exp_Bp
    sigma_Bp_eff = cc_out.sigma_Bp * cc_in.sigma_Bp
    sigma_RphZ_eff = cc_out.sigma_RphZ * cc_in.sigma_RphZ
    sigma_rhothph_eff = cc_out.sigma_rhothph * cc_in.sigma_rhothph

    if ~keyword_set(sigma_Ip) then begin
        sigma_Ip_eff = sigma_RphZ_eff
    endif else begin
        sigma_Ip_eff = product(sigma_Ip)
    endelse

    if ~keyword_set(sigma_B0) then begin
        sigma_B0_eff = sigma_RphZ_eff
    endif else begin
        sigma_B0_eff = product(sigma_B0)
    endelse

    g.r = g.r * l_d_eff
    g.rdim = g.rdim * l_d_eff
    g.rleft = g.rleft * l_d_eff
    g.rbbbs = g.rbbbs * l_d_eff
    g.rlim = g.rlim * l_d_eff
    g.rcentr = g.rcentr * l_d_eff
    g.rmaxis = g.rmaxis * l_d_eff
    g.z = g.z * l_d_eff
    g.zdim = g.zdim * l_d_eff
    g.zmid = g.zmid * l_d_eff
    g.zbbbz = g.zbbbz * l_d_eff
    g.zlim = g.zlim * l_d_eff
    g.zmaxis = g.zmaxis * l_d_eff
    g.nbdry = g.nbdry * l_d_eff
    g.lim = g.lim * l_d_eff

    g.simag = g.simag * sigma_Ip_eff * sigma_Bp_eff * ((2*!dpi)^exp_Bp_eff) * (l_d_eff^2) * l_B_eff
    g.ssimag = g.ssimag * sigma_Ip_eff * sigma_Bp_eff * ((2*!dpi)^exp_Bp_eff) * (l_d_eff^2) * l_B_eff
    g.sibry = g.sibry * sigma_Ip_eff * sigma_Bp_eff * ((2*!dpi)^exp_Bp_eff) * (l_d_eff^2) * l_B_eff
    g.ssibry = g.ssibry * sigma_Ip_eff * sigma_Bp_eff * ((2*!dpi)^exp_Bp_eff) * (l_d_eff^2) * l_B_eff
    g.psi = g.psi * sigma_Ip_eff * sigma_Bp_eff * ((2*!dpi)^exp_Bp_eff) * (l_d_eff^2) * l_B_eff
    g.psirz = g.psirz * sigma_Ip_eff * sigma_Bp_eff * ((2*!dpi)^exp_Bp_eff) * (l_d_eff^2) * l_B_eff
    
    g.bcentr = g.bcentr * l_B_eff * sigma_B0_eff
    g.current = g.current * sigma_Ip_eff * l_d_eff * l_B_eff / (mu0^exp_mu0_eff)
    g.fpol = g.fpol * sigma_B0_eff * l_d_eff * l_B_eff
    g.pres = g.pres * (l_d_eff^2) / (mu0^exp_mu0_eff)
    g.ffprim = g.ffprim * sigma_Ip_eff * sigma_Bp_Eff / ((2*!dpi)^exp_Bp_eff) * l_B_eff
    g.pprime = g.pprime * sigma_Ip_eff * sigma_Bp_Eff / ((2*!dpi)^exp_Bp_eff) * l_B_eff / ((mu0^exp_mu0_eff) * (l_D_eff^2))
    g.qpsi = g.qpsi * sigma_Ip_eff * sigma_B0_eff * sigma_rhothph_eff
END
