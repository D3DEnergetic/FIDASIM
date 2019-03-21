FUNCTION test_beam,beta

    uvw_src = [530.d0+2.d0*cos(beta), 0.d0, 2.d0*sin(beta)]
    uvw_pos = [530.d0, 0.d0, 0.d0]
    uvw_axis = uvw_pos - uvw_src
    uvw_axis = uvw_axis/sqrt(total(uvw_axis*uvw_axis))

    focy = 999999.9d0
    focz = 1000d0
    divy = replicate(8.73d-3,3)
    divz = replicate(2.27d-2,3)
    widy=6d0
    widz=24d0
    naperture = 1
    ashape = [1]
    awidy = [8.85d0]
    awidz = [24d0]
    aoffy = [0.d0]
    aoffz = [0.d0]
    adist = [186.1d0]
    nbi={name:'test_beam',shape:1,data_source:'test_beam.pro',$
         src:uvw_src,axis:uvw_axis,widy:widy,widz:widz,$
         divy:divy,divz:divz,focy:focy,focz:focz, $
         naperture:naperture,ashape:ashape,adist:adist, $
         awidy:awidy,awidz:awidz,aoffy:aoffy,aoffz:aoffz }

    return, nbi
END
