;======================================================================
; $Id: calculate_bfield.pro,v 1.3 2008/06/14 18:51:54 liud Exp $
;======================================================================
;  WWH removed requirement for a eqdsk from Jong's calculate_bfields
; Include sign of current so br & bz go in right direction for 
; standard (r,phi,z) coordinates
; NAME: 
;     CALCULATE_BFIELD
;
; PURPOSE: 
;     Calculate the poloidal, toroidal, vertical and radial magnetic
;     field components
;
; CALLING SEQUENCE: 
;     calculate_bfields,bp,br,bt,bz,a,g
;
; INPUT PARAMETERS: 
;     a:      - structure containing A0 parameters
;     g:      - structure containing G0 parameters
;
;
; OPTIONAL INPUT PARAMETERS: 
;
;     NONE
;
; KEYWORDS: 
;
;     NONE
;
; OUTPUTS: 
;      bp     2-D array containing the poloidal field at refit,zefit
;      br     2-D array containing the radial field at refit,zefit
;      bt     2-D array containing the toroidal field at refit,zefit
;      bz     2-D array containing the vertical field at refit,zefit
;
; COMMON BLOCKS: 
;
;     NONE
;
; SIDE EFFECTS: 
;
;     NONE
;
; RESTRICTIONS:
;     Prior to this using this function, you must have first read in
;     the EFIT data to fill the a and g structures. 
;
; PROCEDURE: 
;     The program uses a simple differentiation using a 3 point
;     Lagrangian interpolation.
;
; CODE TYPE: modeling, analysis
;
; CODE SUBJECT:  edge, transport, equilibrium
;
; EASE OF USE: can be used with existing documentation
;
; OPERATING SYSTEMS:  Unix Of All Flavors
;
; EXTERNAL CALLS:  NONE
;
; RESPONSIBLE PERSON: Ray Jong
;	
; DATE OF LAST MODIFICATION:  09/24/98
;
; MODIFICATION HISTORY:
;
;     Created by Gary D. Porter, LLNL
;     1994.02.16     Michael D. Brown
;                    Optimized for IDL
;     1998.02.17:    Gary D. Porter
;                    Modified to use new EFIT routines and a and g structures.
;                    No longer uses efitcommon.
;-	

PRO calculate_bfield,bp,br,bt,bz,g

    compile_opt defint32,strictarr,strictarrsubs
    
    if n_elements(g.time) eq 0 then return
    ;if a.ishot le 0 or g.time le 0 then return
    if g.time le 0 then return
    mw=g.mw & mh=g.mh
    bp=fltarr(mw,mh) & bt=fltarr(mw,mh) & br=fltarr(mw,mh) & bz=fltarr(mw,mh)
    dpsidx = fltarr(mw,mh)
    dpsidy = fltarr(mw,mh)
    
    ; calculate vertical derivative of psi
    for i = 0,mw-1 do begin
     dpsidy[i,*] = Deriv(g.z[0:mh-1],g.psirz[i,0:mh-1])
    endfor
    
    ; calculate horizontal derivative of psi
    for j = 0,mh-1 do begin
      dpsidx[*,j] = Deriv(g.r[0:mw-1],g.psirz[0:mw-1,j])
    endfor
    
    ; calculate array of Br, Bz, and Bp
    for j = 0,mh-1 do begin
       br[*,j] = dpsidy[0:mw-1,j]/g.r[0:mw-1]
       bz[*,j] = -dpsidx[0:mw-1,j]/g.r[0:mw-1]
    endfor
    bp = sqrt(br*br+bz*bz)
    
    ; WWH get right sign
    if g.cpasma lt 0. then begin
      br=-br & bz=-bz
    end
    
    ; Calculate toroidal field
    ; Original coding was from gfield.for by Peter Politzer,
    ;   translated to IDL by Gary Porter (see BFIELD.PRO).
    ; The code below has be optimized for IDL by Michael D. Brown, 2/16/94
    
    dpsi = (g.ssibry-g.ssimag)/float(mw-1)
    ; first order Bt value.
    for j=0,mh-1 do bt[0:mw-1,j]=g.bcentr*g.rzero/g.r[0:mw-1]  
    k = long((g.psirz - g.ssimag)/dpsi)
    iw=where(k ge 0 and k lt mw-1,n)  ; 1-d indexes where k is a valid index.
    if n gt 0 then begin
      iwr = iw mod mw  ; map matrix 1-d selected indexes to an refit row index.
      bt[iw] = ( g.fpol[k[iw]]+(g.fpol[k[iw]+1]-g.fpol[k[iw]])* $
                 (g.psirz[iw]-(k[iw]*dpsi+g.ssimag))/dpsi ) / g.r[iwr]
    endif
    
    return
end

FUNCTION read_geqdsk,filename,grid,flux=flux
    ;+#read_geqdsk
    ;+Reads an EFIT GEQDSK file
    ;+***
    ;+##Arguments
    ;+    **filename**: GEQDSK file
    ;+
    ;+    **grid**: Interpolation grid
    ;+
    ;+##Keyword Arguments
    ;+    **flux**: Set this keyword to a named variable to recieve the torodial flux mapped onto the interpolation grid
    ;+
    ;+##Return Value
    ;+Electronmagnetic fields structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> fields = read_geqdsk("./g133223.00200",flux=flux)
    ;+```

    equil={err:1}
    ;; Get eqdsk
    g=readg(filename)
    time = double(g.time)
    fluxgrid=double(rho_rz(g,grid.r2d/100.,grid.z2d/100.,/do_linear))

    calculate_bfield,bp,br,bphi,bz1,g

    ;; Get radial electric field on efit's grid from potential
    ;; epoten is on a grid of equally spaced points in psi from g.ssimag to g.ssibry
    dpsi=(g.ssibry-g.ssimag)/(n_elements(g.epoten)-1)
    psi=g.ssimag + dpsi*findgen(n_elements(g.epoten))
    npot=n_elements(g.epoten)
    epot=replicate(g.epoten[npot-1],n_elements(g.r),n_elements(g.z))
    for i=0l,n_elements(g.r)-1 do begin
        for j=0l,n_elements(g.z)-1 do begin
            psi1=g.psirz[i,j]
                        dum=min(abs(psi1-psi),kpsi)
                        if kpsi ne npot-1 then epot[i,j]=spline(psi,g.epoten,[psi1])
              endfor
    endfor

    ; E = - grad(Phi)    EFIT units should be V/m
    er=-(shift(epot,-1,0) - shift(epot,1,0))/(g.r[2]-g.r[0])
    ez1=-(shift(epot,0,-1) - shift(epot,0,1))/(g.z[2]-g.z[0])
  
    ;; Interpolate cylindrical fields onto (r,w) mesh
    b_r=dblarr(grid.nr,grid.nz) & b_t=b_r & b_z=b_r
    e_r=dblarr(grid.nr,grid.nz) & e_t=e_r & e_z=e_r

    for i=0L,grid.nr-1 do for j=0L,grid.nz-1 do begin
            rgrid=(.01*grid.r2d[i,j] - g.r[0])/(g.r[1]-g.r[0]) ; in grid units
            zgrid=(.01*grid.z2d[i,j] - g.z[0])/(g.z[1]-g.z[0])    ; WWH 3/31/07
            b_r[i,j]  =interpolate(br,[rgrid],[zgrid])
            e_r[i,j]  =interpolate(er,[rgrid],[zgrid])
            b_t[i,j]  =interpolate(bphi,[rgrid],[zgrid])
            e_z[i,j]  =interpolate(ez1,[rgrid],[zgrid])
            b_z[i,j]  =interpolate(bz1,[rgrid],[zgrid])
    endfor
  
    flux = fluxgrid
    mask = replicate(1,grid.nr,grid.nz)

    equil={time:time,data_source:filename, mask:mask, $
           br:b_r,bt:b_t,bz:b_z,er:e_r,et:e_t,ez:e_z}
    GET_OUT:
    return,equil
END
