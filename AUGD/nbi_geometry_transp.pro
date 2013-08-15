pro rotate_uvw,uvw,Arot,Brot,Crot,updown,xyz
  A2rot=reform(Arot[0,*,*])
  B2rot=reform(Brot[0,*,*])
  C2rot=reform(Crot[0,*,*])
  ;;rotate uvz by alpa alround z axis
  if updown lt 0 then qrz=MATRIX_MULTIPLY(A2rot,uvw)
  if updown ge 0 then qrz=MATRIX_MULTIPLY(B2rot,uvw)
  ;; rotate qrz_ray by phi_box on xyz_coordinates
  xyz=MATRIX_MULTIPLY(C2rot,qrz)
end

PRO nbi_geometry_transp,nbgeom, doplt=doplt,do_oplot_tor=do_oplot_tor $
                        ,rotate=rotate,src=src,cm=cm
;===============================================
; Beam geometry settings from AUGD homepage
;===============================================
  if not keyword_set(shot_in) then shot=27237 else shot=shot_in
;;R(P):: Distance horizontal beam crossing to - torus axis [cm]
  R0=[284.20d0,284.20d0,284.20d0,284.20d0,329.63d0,329.63d0,329.63d0,329.63d0]
;;PHI: angle between R and box [rad]
  PHI=[15.d0,15.d0,15.d0,15.d0,18.9d0,18.9d0,18.9d0,18.9d0]/!radeg
;;THETA:   angle towards P (horizontal beam crossing) [rad]
  THETA=[33.75d0,33.75d0,33.75d0,33.75d0,29.d0,29.d0,29.d0,29.d0]/!radeg
  THETA[4:7]=  THETA[4:7] + !pi   ;; for NBI box2
;;ALPHA: horizontal angle between Box-axis and source [rad]
  alpha=4.1357d0/!radeg
  ALPHA=alpha*[1.d0,-1.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0]
;;BETA: vertical angle between box-axis and source [rad]
 beta  = 4.8991d0/!radeg
 beta6 = 6.65d0/!radeg 
 beta7 = 6.65d0/!radeg 
 BETA  = [-beta,-beta,beta,beta,-beta,-beta6,beta7,beta]
 
;;Ion source half width [cm]
  bmwidra=11.d0 
  bmwidza=25.d0 
;;Divergence of beams [rad]
  DIVZA=replicate(0.0111,3,8)
  DIVRA=replicate(0.0111,3,8)

;;distance source P0
  d_src_P0=650.0d0
;;Focal lengths [cm]
  foclra=replicate(650.0d0,8) ;horizontal
  foclza=replicate(850.0d0,8) ;vertical
;;Half horizontal displacement of sources (e.g. #1, #3) [cm]
  src_hw=47.d0 
  src_hw=src_hw*[1.d0,-1.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0] 
;;source elevation w.r.t. midplane [cm]
  src_hv=[60.d0,60.d0,-60.d0,-60.d0,60.d0,70.d0,-70.d0,-60.d0]
;;Width of the aperture in point P [cm]
  rapedga=20.d0
  xzpedga=30.d0 
;;Phi_box: angle of the box-axi to the x-axis [rad]
  phi_box = THETA-PHI
;;x/ypos_hcrossing: position of horizontal beam crossing [cm]
  xpos_P0 = R0 * cos(THETA)
  ypos_P0 = R0 * sin(THETA)
  zpos_P0 = ypos_P0*0.d0
  xyz_P0=dblarr(8,3)
  for k=0,7 do begin
     xyz_P0[k,*]=[xpos_P0[k],ypos_P0[k], zpos_P0[k]]
  endfor
;;x/zpos_vcrossing: position of vertical beam crossing [cm]
  delta_u=[-50.d0,-50.d0,-50.d0,-50.d0,50.d0,50.d0,50.d0,50.d0]
  delta_v=[-3.62d0,+3.62d0,+3.62d0,-3.62d0,0,0,0,0]
  xpos_P1= xpos_P0+delta_u*cos(-phi_box)+sin(-phi_box)*delta_v
  ypos_P1= ypos_P0-delta_u*sin(-phi_box)+cos(-phi_box)*delta_v
  zpos_P1= 0.d0    





;;x/y/zpos_source: position of source in xy coordinates
  xpos_source=xpos_P0 +  d_src_P0*cos(phi_box+ALPHA)
  ypos_source=ypos_P0 +  d_src_P0*sin(phi_box+ALPHA)
  zpos_source=src_hv 
  xyz_src=dblarr(8,3)
  for k=0,7 do begin
     xyz_src[k,*]=[xpos_source[k],ypos_source[k], zpos_source[k]]
  endfor
  xyz_pos=dblarr(8,3)
  xyz_pos[*,0]=xpos_P0
  xyz_pos[*,1]=ypos_P0
  xyz_pos[*,2]=zpos_P0
  xyz_plasma=dblarr(8,3)
  xyz_vec=dblarr(8,3) 
  for ii=0,7 do begin
     xyz_vec[ii,*]=xyz_pos[ii,*]-xyz_src[ii,*]
     xyz_vec[ii,*]=xyz_vec[ii,*] $
                   /sqrt(xyz_vec[ii,0]^2+xyz_vec[ii,1]^2+xyz_vec[ii,2]^2)
     xyz_plasma[ii,*]=xyz_src[ii,*]+xyz_vec[ii,*]*1000.d0
  endfor


  gamma=-!pi/8.*3.
  Rrot=dblarr(3,3)
  Rrot[0,0]= cos(GAMMA)  & Rrot[0,1]=-sin(GAMMA)  & Rrot[0,2]= 0.
  Rrot[1,0]= sin(GAMMA)  & Rrot[1,1]= cos(GAMMA)  & Rrot[1,2]= 0.
  Rrot[2,0]= 0.          & Rrot[2,1]= 0.          & Rrot[2,2]= 1. 
  for ii=0,7 do begin
     xyz_src[ii,*]=MATRIX_MULTIPLY(Rrot, reform(xyz_src[ii,*]))
     xyz_pos[ii,*]=MATRIX_MULTIPLY(Rrot, reform(xyz_pos[ii,*]))
     xyz_plasma[ii,*]=MATRIX_MULTIPLY(Rrot, reform(xyz_plasma[ii,*]))
  endfor
  phi_box=phi_box+replicate(gamma,8)
  

  ;; If another coordinate system is used (e.g. by rmm)
  if keyword_set(rotate) then begin
     Rrot=dblarr(3,3)
     Rrot[0,0]= cos(ROTATE)  & Rrot[0,1]=-sin(ROTATE)  & Rrot[0,2]= 0.
     Rrot[1,0]= sin(ROTATE)  & Rrot[1,1]= cos(ROTATE)  & Rrot[1,2]= 0.
     Rrot[2,0]= 0.          & Rrot[2,1]= 0.          & Rrot[2,2]= 1. 
     for ii=0,7 do begin
        xyz_src[ii,*]=MATRIX_MULTIPLY(Rrot, reform(xyz_src[ii,*]))
        xyz_pos[ii,*]=MATRIX_MULTIPLY(Rrot, reform(xyz_pos[ii,*]))
        xyz_plasma[ii,*]=MATRIX_MULTIPLY(Rrot, reform(xyz_plasma[ii,*]))
     endfor
     phi_box=phi_box+replicate(rotate,8)
  endif




  zero=replicate(0.d0,8)
  one=replicate(1.d0,8)
;;transformation matrix to rotate on NBI box axis by alpha
  alpha=phi_box+alpha
  Arot=dblarr(8,3,3)
  Arot[*,0,0]= cos(BETA)   & Arot[*,0,1]= zero    & Arot[*,0,2]= sin(BETA)
  Arot[*,1,0]= zero        & Arot[*,1,1]= one    & Arot[*,1,2]= zero 
  Arot[*,2,0]=-sin(BETA)   & Arot[*,2,1]= zero    & Arot[*,2,2]= cos(BETA) 
  ;; transformation matrix to rotate in vertical direction by beta
  Brot=dblarr(8,3,3)
  Brot[*,0,0]= cos(BETA)   & Brot[*,0,1]= zero    & Brot[*,0,2]= sin(BETA)
  Brot[*,1,0]= zero        & Brot[*,1,1]= one    & Brot[*,1,2]= zero 
  Brot[*,2,0]=-sin(BETA)   & Brot[*,2,1]= zero    & Brot[*,2,2]= cos(BETA) 
  ;; transformation matrix to rotate towards the xy-coordinate system
  Crot=dblarr(8,3,3)
  Crot[*,0,0]= cos(alpha) & Crot[*,0,1]=-sin(alpha) & Crot[*,0,2]= zero
  Crot[*,1,0]= sin(alpha) & Crot[*,1,1]= cos(alpha) & Crot[*,1,2]= zero 
  Crot[*,2,0]= zero         & Crot[*,2,1]= zero         & Crot[*,2,2]= one




;=========================
; store data in structure
;=========================
  nbgeom={  Arot:   Arot $
             , Brot:   Brot $
             , Crot:   Crot $
             ;; parameters
             , focy: double(foclra)     , focz: double(foclza) $
             , divy: double(divra)      , divz: double(divza) $
             , bmwidra:double(bmwidra)  , bmwidza:double(bmwidza)  $
             , xyz_src: xyz_src, xyz_pos: xyz_pos  }
    

;==========
; Graphics
;==========
  if keyword_set(doplt) then begin 
     loadct,39   ;; green=150;blue=50;yellow=200;red=254;black=0;white=255
     device, decomposed=0
     !P.background=255 & !P.color=0 &!P.font=-1 
     !P.charsize=1.6 & !p.multi=[0,0,0]
     xran=[-300.,900]
     yran=[-300.,500]
     zran=[-100.,100]
     if not(keyword_set(ps)) then  window,1,colors=15,xsize=500,ysize=500
     plot,[0.],/nodata,xrange=xran,yrange=yran,xtit='X [cm]',ytit='Y [cm]' $
          ,/isotropic
      rmaj=165. &  rmin=50. &  n_th=100 
      theta_asdex=findgen(n_th+1)*2.*!PI/(n_th)
      r_in =fltarr(n_th+1)+rmaj-rmin
      r_out=fltarr(n_th+1)+rmaj+rmin
      oplot,[0.],[0.],psym=1,thick=2
      oplot, r_in , theta_asdex,/polar
      oplot,rmaj+fltarr(n_th+1),theta_asdex,/polar,linestyle=2
      oplot, r_out, theta_asdex,/polar  
      for ii=0,7 do begin
         oplot,[xyz_pos[ii,0]],[xyz_pos[ii,1]] ,psym=1,color=200
         oplot,[xyz_pos[ii,0],xyz_src[ii,0]] $
               ,[xyz_pos[ii,1],xyz_src[ii,1]],color=254

         oplot,[xyz_src[ii,0],xyz_plasma[ii,0]] $
               ,[xyz_src[ii,1],xyz_plasma[ii,1]],color=254

         ;;plot rays      
         uvw_ray=[-1.d0,0.d0,0.d0]*1000.
         updown=1
         rotate_uvw,uvw_ray,Arot[ii,*,*],Brot[ii,*,*],Crot[ii,*,*],updown $
                    ,xyz_ray
         oplot,[xyz_src[ii,0],xyz_src[ii,0]+xyz_ray[0]] $
               ,[xyz_src[ii,1],xyz_src[ii,1]+xyz_ray[1]],thick=2,color=50
      endfor


  if not(keyword_set(ps)) then  window,2,colors=15,xsize=1200,ysize=400
 plot,[0.],/nodata,xrange=[00,900],yrange=zran,xtit='X [cm]',ytit='Z [cm]' $
      ,/isotropic
 oplot,[0.],[0.],psym=1,thick=2
 x_torus=[-rmaj-rmin,-rmaj,-rmaj+rmin,rmaj-rmin,rmaj,rmaj+rmin]
 y_torus=[-1.6*rmin,1.6*rmin]
 for jj=0,5 DO BEGIN
    sty=0
    if (ABS(x_torus(jj)) eq rmaj) then sty=3
         oplot,[x_torus(jj),x_torus(jj)],[y_torus(0),y_torus(1)],linestyle=sty
      endfor
 for ii=0,3 do begin
    oplot,[xyz_pos[ii,0]],[xyz_pos[ii,2]] ,psym=1
    oplot,[xyz_pos[ii,0],xyz_src[ii,0]] $
          ,[xyz_pos[ii,2],xyz_src[ii,2]],color=254
    oplot,[xyz_pos[ii,0],xyz_plasma[ii,0]] $
          ,[xyz_pos[ii,2],xyz_plasma[ii,2]],color=254    
    oplot,[xyz_src[ii,0]], [xyz_src[ii,2]],color=1,psym=4
    
    ;; plot rays
    uvw_ray=[-1.d0,0.d0,0.d0]*1000.
    updown=-1
    rotate_uvw,uvw_ray,Arot[ii,*,*],Brot[ii,*,*],Crot[ii,*,*],updown,xyz_ray
    oplot, [xyz_src[ii,0],xyz_src[ii,0]+xyz_ray[0]],$
           [xyz_src[ii,2],xyz_src[ii,2]+xyz_ray[2]],color=50 

 endfor
   endif


  if keyword_set(do_oplot_tor) then begin
     for ii=0,7 do begin
        if keyword_set(src) then if ii ne src then continue
                ;;plot rays      
         uvw_ray=[-1.d0,0.d0,0.d0]*1000.
         updown=1
         rotate_uvw,uvw_ray,Arot[ii,*,*],Brot[ii,*,*],Crot[ii,*,*],updown $
                    ,xyz_ray
         if keyword_set(cm) then begin
            xyz_ray=xyz_ray*100.
            xyz_src=xyz_src*100.
         endif
         oplot,[xyz_src[ii,0],xyz_src[ii,0]+xyz_ray[0]]/100. $
               ,[xyz_src[ii,1],xyz_src[ii,1]+xyz_ray[1]]/100.,thick=2,color=0
        
     endfor
  endif
END
