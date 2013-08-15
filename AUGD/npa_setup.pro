pro npa_setup,det
      nchan=1
      xlos= [136.]
      ylos= [-25.]
      zlos=  [10.]
    ;; position at NPA
      uhead = 0.95*100. ;[cm]
      vhead = 3.32*100. ;[cm]
      ;; turn head position by 360/16*3 deg
      phihead=atan(vhead/uhead)
      rhead=sqrt(vhead^2+uhead^2)
      phihead=phihead-(360.d0/16.d0*3.d0)/!radeg
      xhead=[rhead*cos(phihead)]
      yhead=[rhead*sin(phihead)]
      zhead = [10.]                   ;[cm]
      
      ;; THE headsize and opening angle determine the flux of neutrals
      ;; on the npa detector!
      headsize = [0.6] ;; cm
      opening_angle = [1./180.*!pi]

      det = set_det( NCHAN=nchan, $
                     XLOS=xlos, YLOS=ylos, ZLOS=zlos, $
                     XHEAD=xhead, YHEAD=yhead, ZHEAD=zhead, $
                     HEADSIZE=headsize, OPENING_ANGLE=opening_angle)
end
