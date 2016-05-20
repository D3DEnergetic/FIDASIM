FUNCTION tb_zyx, a, b, g
    ;; Returns Tait-Bryan z-y'-x" active rotation matrix
    ;; a: (alpha) rotation about z
    ;; b: (beta)  rotation about y'
    ;; g: (gamma) rotation about x"
    sa = sin(a) & ca = cos(a)
    sb = sin(b) & cb = cos(b)
    sg = sin(g) & cg = cos(g)

    R = dblarr(3,3)
    R[0,0] = ca*cb & R[1,0] = ca*sb*sg - cg*sa & R[2,0] = sa*sg + ca*cg*sb
    R[0,1] = cb*sa & R[1,1] = ca*cg + sa*sb*sg & R[2,1] = cg*sa*sb - ca*sg
    R[0,2] = -sb   & R[1,2] = cb*sg            & R[2,2] = cb*cg

    return, R
END

FUNCTION line_to_rotmat,r0,rf

    dis = sqrt(total((rf - r0)^2.0))
    alpha = atan(rf[1] - r0[1],rf[0] - r0[0])
    beta = asin((r0[2] - rf[2])/dis)
    return, tb_zyx(alpha,beta,0.d0)

END

FUNCTION test_npa
    
    ;; Chords
    nchan = 3L
    ulens = dblarr(nchan)
    vlens = [-170.d0,-170.d0,-170.d0]
    wlens = replicate(100.d0,nchan)
    lens = transpose([[ulens],[vlens],[wlens]])
    ulos = dblarr(nchan)
    vlos = [-200.d0,-170.d0,-140.d0]
    wlos = dblarr(nchan)
    radius = sqrt(ulos^2.d0 + vlos^2.d0)
    id = ["c1","c2","c3"]

    a_cent  = dblarr(3,nchan)
    a_redge = dblarr(3,nchan)
    a_tedge = dblarr(3,nchan)
    d_cent  = dblarr(3,nchan)
    d_redge = dblarr(3,nchan)
    d_tedge = dblarr(3,nchan)

    ac = [0.d0, 0.d0, 0.d0]
    ar = [0.d0, 3.d0, 0.d0]
    at = [0.d0, 0.d0, 3.d0]

    dc = [-50.d0, 0.d0, 0.d0]
    dr = [-50.d0, 3.d0, 0.d0]
    dt = [-50.d0, 0.d0, 3.d0]
    
    for i=0,2 do begin
        r0 = [ulens[i],vlens[i],wlens[i]]
        rf = [ulos[i],vlos[i],wlos[i]]
        R = line_to_rotmat(r0,rf)
 
        a_cent[*,i] = transpose(R##transpose(ac)) + r0
        a_redge[*,i] = transpose(R##transpose(ar)) + r0
        a_tedge[*,i] = transpose(R##transpose(at)) + r0
        
        d_cent[*,i] = transpose(R##transpose(dc)) + r0
        d_redge[*,i] = transpose(R##transpose(dr)) + r0
        d_tedge[*,i] = transpose(R##transpose(dt)) + r0
    endfor

    npa_chords = {nchan:3L,system:"NPA",data_source:"test_npa.pro", id:id, $
                  a_shape:replicate(2,nchan),d_shape:replicate(2,nchan), $
                  a_cent:a_cent,a_redge:a_redge,a_tedge:a_tedge, $
                  d_cent:d_cent,d_redge:d_redge,d_tedge:d_tedge, radius:radius}

    return, npa_chords
END
