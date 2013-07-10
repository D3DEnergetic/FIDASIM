pro load_bfield,inputs,coords,expnam,dianam,b,rhoPF,rhoTF
  ;; Magnetic field vectors and rho from Cliste Equilibrium at ASDEX Upgrade
  ;; INPUTS: inputs.shot
  ;;         inputs.time
  ;;         coords.ng
  ;;         coords.r_grid
  ;;         coords.zc
  ;;         dianam (name of Equilibrium)
  ;; OUTPUT: MAGNETIC field vector in xyz coordinates
  ;;         RHO poloidal
  ;;         RHO toroidal
  nSHOT=long(inputs.shot)
  nEDIT=0L
  tSHOT=float(inputs.time)
  IERR=0L 
  rin=float(coords.r_grid/100.)
  zin=float(coords.zc/100.)
  lin=long(coords.ng)

  ;; LOAD magnetic field vectors
  Br=fltarr(lin)
  Bz=fltarr(lin)
  Bt=fltarr(lin)
  fPF=fltarr(lin)
  fTF=fltarr(lin)
  s = call_external( !libkk, 'kkidl', 'kkrzBrzt' $
                     , iERR  ,expnam,dianam,nSHOT,nEDIT ,tSHOT $
                     , rin, zin, lin $
                     , Br, Bz, Bt, fPF, fTF)      
  b=dblarr(3,coords.ng)
  b[0,*] =- cos(!pi*0.5d0-coords.phi_grid)*Bt + cos(coords.phi_grid)*Br
  b[1,*] =  sin(!pi*0.5d0-coords.phi_grid)*Bt + sin(coords.phi_grid)*Br 
  b[2,*] =  Bz  

  ;; LOAD grid for flux coordinates (rho)
  fPF=fltarr(lin)
  fTF=fltarr(lin)
  rhoPF=fltarr(lin)
  rhoTF=fltarr(lin)
  s = CALL_EXTERNAL(!libkk,'kkidl','kkrzPTFn' $
                    , iERR, expnam,dianam,nSHOT, nEDIT, tSHOT $
                    , rin, zin, lin, fPF, rhoPF, fTF, rhoTF)
  
  ;; extent rhoTF by rhoPF outside the separatri
  index=where(rhoPF gt 1.)
  rhoTF[index]=rhoPF[index]
end
