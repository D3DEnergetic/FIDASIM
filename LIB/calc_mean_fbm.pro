PRO calc_mean_fbm,grid,fbm,weights,los,neutrals,mean_fbm2,elevel=elevel

	if not keyword_set(elevel) then elevel=indgen(n_elements(neutrals.fdens[0,0,0,*]))

	nchan=los.nchan
    nen_transp=fbm.fbm_nenergy
    npitch_transp=fbm.FBM_npitch
    dE_transp=fbm.fbm_energy[1]-fbm.fbm_energy[0]
    dpitch_transp=abs(fbm.fbm_pitch[1]-fbm.fbm_pitch[0])

    r2d=fbm.fbm_r2d
    z2d=fbm.fbm_z2d
    full_fbm=fbm.fbm
    mean_fbm= fltarr(nen_transp,npitch_transp,nchan)
    mean_fbm2=fltarr(n_elements(weights.energy),n_elements(weights.pitch),nchan)

    for ichan=0, nchan-1 do begin
        ;; ------------------------------------------
        ;;------------ CALCUATE mean value of fbm----
        ;; ------------------------------------------
        rad=0.d0
        for i=0,grid.nx-1 do begin
           for j=0,grid.ny-1 do begin
              for k=0,grid.nz-1 do begin
                 los_wght=los.los_wght[i,j,k,ichan]
                 if los_wght gt 0. then begin
                    ;; determine mean values like the halo density along LOS
                    wght=(  total(neutrals.fdens[i,j,k,elevel] + $
                            neutrals.hdens[i,j,k,elevel] + $
                            neutrals.tdens[i,j,k,elevel] + $
                            neutrals.halodens[i,j,k,elevel])) * los_wght
                    rad=rad+wght
					rrc=grid.r_grid[i,j,k]
					zzc=grid.w_grid[i,j,k]

                    dr=2.       ;[cm]
                    dz=2.       ;[cm]

                    dummy=min((r2d-rrc+dr)^2.+(z2d-zzc+dz)^2.,fbm_index1)
                    dummy=min((r2d-rrc-dr)^2.+(z2d-zzc+dz)^2.,fbm_index2)
                    dummy=min((r2d-rrc+dr)^2.+(z2d-zzc-dz)^2.,fbm_index3)
                    dummy=min((r2d-rrc-dr)^2.+(z2d-zzc-dz)^2.,fbm_index4)
                    dummy=min((r2d-rrc)^2.+(z2d-zzc)^2.,fbm_index5)
                    mean_fbm[*,*,ichan] = mean_fbm[*,*,ichan] + ( $
                                       full_fbm[*,*,fbm_index1] + $
                                       full_fbm[*,*,fbm_index2] + $
                                       full_fbm[*,*,fbm_index3] + $
                                       full_fbm[*,*,fbm_index4] + $
                                       full_fbm[*,*,fbm_index5])/5.*wght 
                 endif
              endfor
           endfor
        endfor
        mean_fbm[*,*,ichan]=mean_fbm[*,*,ichan]/rad
        ;;------------------------------------------------------------
        ;; map FBM on the energy and pitch grid of the weight function
        ;;------------------------------------------------------------
        for ie=0,n_elements(weights.energy)-1 do begin
           dummy=min(abs(fbm.fbm_energy-weights.energy[ie]),eindex)
           for ip=0,n_elements(weights.pitch)-1 do begin
              dummy=min(abs(fbm.fbm_pitch-weights.pitch[ip]),pindex)
              mean_fbm2[ie,ip,ichan]=mean_fbm[eindex,pindex,ichan]
              ;;[fast_ion/cm^3/dP/dE]
           endfor
        endfor
  	endfor ;; loop over channels	
END
