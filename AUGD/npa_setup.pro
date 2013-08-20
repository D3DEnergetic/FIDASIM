FUNCTION set_det,nchan=nchan,xlos=xlos,ylos=ylos,zlos=zlos,xhead=xhead,yhead=yhead, $
                 zhead=zhead,headsize=headsize,opening_angle=opening_angle, $
                 sigma_pi=sigma_pi,nx=nx,ny=ny,xx=xx,yy=yy,xyzpix=xyzpix

  IF NOT KEYWORD_SET(nchan) THEN nchan = 1
  IF NOT KEYWORD_SET(xlos)  THEN xlos = [0.]
  IF NOT KEYWORD_SET(ylos)  THEN ylos = [0.] 
  IF NOT KEYWORD_SET(zlos)  THEN zlos = [0.]
  IF NOT KEYWORD_SET(xhead) THEN xhead = 0.
  IF NOT KEYWORD_SET(yhead) THEN yhead = 0.
  IF NOT KEYWORD_SET(zhead) THEN zhead = 0.
  IF NOT KEYWORD_SET(opening_angle) THEN opening_angle = dblarr(nchan)
  IF NOT KEYWORD_SET(headsize)      THEN headsize = dblarr(nchan)
  IF NOT KEYWORD_SET(sigma_pi)      THEN sigma_pi = dblarr(nchan)+1.0

  det = { nchan: nchan, $
          xlos: xlos, ylos: ylos, zlos: zlos, $
          xhead: xhead, yhead: yhead, zhead: zhead, $
          npa: {opening_angle: opening_angle, headsize: headsize }, $
          sigma_pi: sigma_pi }
  RETURN,det
END

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
