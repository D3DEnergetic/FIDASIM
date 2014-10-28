;======================================================================
; $Id: calculate_bfield.pro,v 1.3 2008/06/14 18:51:54 liud Exp $
;======================================================================
;  WWH removed requirement for a eqdsk from Jong's calculate_bfields
; Include sign of current so br & bz go in right direction for 
; standard (r,phi,z) coordinates
;+ 
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

pro calculate_bfield,bp,br,bt,bz,g

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
