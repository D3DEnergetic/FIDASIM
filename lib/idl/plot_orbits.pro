PRO plot_orbits, orb, g
    ;+#plot_orbits
    ;+Plots charged fusion product orbits
    ;+***
    ;+##Arguments
    ;+    **orb**: Orbits structure (table)
    ;+
    ;+    **g**: GEQDSK file
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> restore, 'mast.idl'
    ;+IDL> g = 'g99999K26'
    ;+IDL> plot_orbits, orb, g
    ;+```

    device, decomposed=0

    g=readg(g)
    finewall,g,rwall,zwall
    sightline = orb.sightline
    daomega = orb.daomega
    nactual = orb.nactual
    nch = 4

    !p.multi=[0,2,0]
    !p.background=255 ;white background
    ; (R,z) elevation

    loadct, 39
    plot,g.lim(0,*),g.lim(1,*),color=0, xtitle='R [m]',ytitle='Z [m]',charsize=2
    oplot,g.bdry(0,0:g.nbdry-1),g.bdry(1,0:g.nbdry-1),color=40

    colors=0 + 80*indgen(4)
    nrays = orb[0].nrays
    for ich=0,nch-1 do for iray=0,nrays-1 do if daomega[iray,ich] gt 0. then begin
        nact=nactual[iray,ich]-1
        oplot,sightline[3,0:nact,iray,ich],sightline[5,0:nact,iray,ich],psym=3,color=colors[ich]
    end

    ; (R,phi) plan
    theta=2.*!pi*findgen(31)/30.

    plot,max(rwall)*cos(theta),max(rwall)*sin(theta),xrange=[-2.5,2.5],yrange=[-4,4],xstyle=1,$
         xtitle='X [m]',ytitle='Y [m]',color=0,charsize=2
    oplot,min(rwall)*cos(theta),min(rwall)*sin(theta),color=0

    for ich=0,nch-1 do for iray=0,nrays-1 do if daomega[iray,ich] gt 0. then begin
         nact=nactual[iray,ich]-1
         oplot,sightline[3,0:nact,iray,ich]*cos(sightline[4,0:nact,iray,ich]), $
               sightline[3,0:nact,iray,ich]*sin(sightline[4,0:nact,iray,ich]),color=colors[ich]
    end


END
