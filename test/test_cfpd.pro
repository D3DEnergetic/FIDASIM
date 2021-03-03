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

FUNCTION test_cfpd, tname, gname
    ;; Chords
    nchan = 4L

    ;Positions
    ZDist=fltarr(4) ; [m]
        ZDist(1-1)= -1.30
        ZDist(2-1)= -1.30
        ZDist(3-1)= -1.30
        ZDist(4-1)= -1.30

    RDist=fltarr(4) ;[m]
        RDist(1-1)= 1.1
        RDist(2-1)= 1.2
        RDist(3-1)= 1.3
        RDist(4-1)= 1.4

    PHDangle=fltarr(4) ;[rad]
        PHDangle(1-1)= 0.d0
        PHDangle(2-1)= 0.d0
        PHDangle(3-1)= 0.d0
        PHDangle(4-1)= 0.d0

    ulens = RDist * cos(PHDangle)
    vlens = RDist * sin(PHDangle)
    wlens = ZDist

    lens = transpose([[ulens],[vlens],[wlens]])

    ;Orientations
    vr = -10 ; radial direction
    vt = +40 ; toroidal direction
    vz = -10 ; vertical direction
    vr /= sqrt(vr^2+vt^2+vz^2)
    vt /= sqrt(vr^2+vt^2+vz^2)
    vz /= sqrt(vr^2+vt^2+vz^2)

    vx = cos(PHDangle)*vr - sin(PHDangle)*vt
    vy = sin(PHDangle)*vr + cos(PHDangle)*vt

    ulos = ulens + vx*(-wlens/vz)
    vlos = vlens + vy*(-wlens/vz)
    wlos = wlens + vz*(-wlens/vz)

    ;Make structure
    radius = sqrt(ulos^2.d0 + vlos^2.d0)
    id = ["c1","c2","c3","c4"]

    a_cent  = dblarr(3,nchan)
    a_redge = dblarr(3,nchan)
    a_tedge = dblarr(3,nchan)
    d_cent  = dblarr(3,nchan)
    d_redge = dblarr(3,nchan)
    d_tedge = dblarr(3,nchan)

    ac = [0.d0, 0.d0, 0.d0]
    ar = [0.d0, 2.5d0, 0.d0]
    at = [0.d0, 0.d0, 2.5d0]

    dc = [-35.7d0, 0.d0, 0.d0]
    dr = [-35.7d0, 2.5d0, 0.d0]
    dt = [-35.7d0, 0.d0, 2.5d0]

    for i=0,3 do begin
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

    cfpd_chords = {nchan:nchan,system:"CFPD",data_source:"test_cfpd.pro", id:id, $
                  a_shape:replicate(2,nchan),d_shape:replicate(2,nchan), $
                  a_cent:a_cent,a_redge:a_redge,a_tedge:a_tedge, $
                  d_cent:d_cent,d_redge:d_redge,d_tedge:d_tedge, radius:radius}

    v=fltarr(3,4)
        v(0,*) = vr
        v(1,*) = vt
        v(2,*) = vz

    D = 0.04      ;                       ! detector-collimator spacing (m)

    RC=fltarr(4)
        RC(1-1) = 1.288e-3  ;                     ! outer collimator radius (m)
        RC(2-1) = 1.294e-3
        RC(3-1) = 1.318e-3
        RC(4-1) = 1.343e-3

    geometry = {rdist:RDist, zdist:ZDist, v:v, d:D, rc:RC}
    earray = 2730. + 150.*findgen(6)
    cfpd_table = proton_table(tname, gname, earray=earray, geometry=geometry, step=0.10, $
                              nsteps=300)
    cfpd = create_struct(cfpd_chords, cfpd_table)

    return, cfpd
END
