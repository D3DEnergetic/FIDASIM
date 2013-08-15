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
