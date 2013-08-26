pro load_bfield,inputs,coords,bx,by,bz,rhoPF,rhoTF
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

  expnam=inputs.equil_exp
  dianam=inputs.equil_diag
  nSHOT=long(inputs.shot)
  nEDIT=long(inputs.equil_ed)
  tSHOT=float(inputs.time)
  IERR=0L 

  lin=long(coords.ng)
  rin=reform(float(coords.rrc_grid/100.),lin)
  zin=reform(float(coords.zzc_grid/100.),lin)

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
  bx =- cos(!pi*0.5d0-coords.phi_grid)*Bt + cos(coords.phi_grid)*Br
  by =  sin(!pi*0.5d0-coords.phi_grid)*Bt + sin(coords.phi_grid)*Br 
  bz =  Bz  
  bx=reform(bx,coords.nx,coords.ny,coords.nz)
  by=reform(by,coords.nx,coords.ny,coords.nz)
  bz=reform(bz,coords.nx,coords.ny,coords.nz)

  ;; LOAD grid for flux coordinates (rho)
  fPF=fltarr(lin)
  fTF=fltarr(lin)
  rhoPF=fltarr(lin)
  rhoTF=fltarr(lin)
  s = CALL_EXTERNAL(!libkk,'kkidl','kkrzPTFn' $
                    , iERR, expnam,dianam,nSHOT, nEDIT, tSHOT $
                    , rin, zin, lin, fPF, rhoPF, fTF, rhoTF)
  ;; reform the arrays !! 
  rhoPF=reform(rhoPF,coords.nx,coords.ny,coords.nz)
  rhoTF=reform(rhoTF,coords.nx,coords.ny,coords.nz) 
  ;; extent rhoTF by rhoPF outside the separatri
  index=where(rhoPF gt 1.)
  rhoTF[index]=rhoPF[index]
end
