;; -----------------------------------------------------------------------------
;;TABLES FOR COLLISIONS BETWEEN PROTONS(H,D) and PROTONS,ELECTRON and IMPURITIES
;; -----------------------------------------------------------------------------
;; 04.2012 by B.Geiger, R.Dux
;; Tables taken from Janev 2004 and ADAS


;; --------------------------------
;; COLLISIONS WITH PROTONS:
;; --------------------------------
;; CHARGE EXCHANGE:
FUNCTION picx_janev,erel,n
  ;; H(+)+H(n)->H(n'l')+H(+) 
  ;; Eq 44 for charge exchange cross section
  ;; in the manual "Collision processes in low-temperature hydrogen plasmas"
  ;; by Janev et al. 
  ;; http://www.eirene.de/report_4105.pdf
  ;; INPUT
  ;; 	erel	protron energy (keV/amu)
  ;;       n       energy level
  ;; OUTPUT        sigma in cm^2
  ;;n=      1           2           3            4
  a1=[ 3.2345,     9.2750e-1,  3.7271e-1,   2.1336e-1 ]
  a2=[ 2.3588e+2,  6.5040e+3,  2.7645e+6,   1.0000e+10]
  a3=[ 2.3713,     2.0699e+1,  1.4857e+3,   1.3426e+6 ]
  a4=[ 3.8371e-2,  1.3405e-2,  1.5720e-3,   1.8184e-3 ]
  a5=[ 3.8068e-6,  3.0842e-6,  3.0842e-6,   3.0842e-6 ]
  a6=[ 1.1832e-10, 1.1832e-10, 1.1832e-10,  1.1832e-10]
  E_=double(n)^2.*Erel
  if n le 4 then i=n-1 else i=3
  sigma=1.e-16*a1[i]*double(n)^4*alog(a2[i]/E_+a3[i])/$
        (1.+a4[i]*E_+a5[i]*E_^3.5+a6[i]*E_^5.4)  
  return, sigma
end
FUNCTION sigma1_adas,erel,mmax
  n=1.
  e=erel*1.e3  ;; conversion from [ev] to [keV]
  ;; set values lower than 1eV to 1ev
  index=where(e lt 1.,nind) 
  if nind gt 0 then e[index]=1.

  nel=n_elements(erel)
  sigma=replicate(0.,mmax,nel)
  ;; m=1
  A=dblarr(20)
  A( 1) =  -3.496092687D+02  
  A( 2) =   4.724931484D+02
  A( 3) =  -2.720493064D+02
  A( 4) =   8.158564625D+01
  A( 5) =  -1.339790721D+01 
  A( 6) =   1.138706949D+00
  A( 7) =  -3.914774156D-02
  A=A[1:7]
  index=where(e ge 1.e3,nind)
  if nind gt 0 then begin
     sigma[0,index]=10^(poly(alog10(e[index]),A))
  endif
  ;; fudge the cross sections for e < 1keV
  index=where(e lt 1.e3,nind)
  if nind gt 0 then begin
     sigma[0,index]=(e[index]*1.e-3)^(-0.2)*10^(poly(alog10(1.e3),A))
  endif
  ;;m=2
  A=dblarr(20)
  A( 1) =  -4.036239511D+03      
  A( 2) =   6.941235312D+03
  A( 3) =  -5.186974866D+03  
  A( 4) =   2.194885201D+03
  A( 5) =  -5.765960509D+02
  A( 6) =   9.653534186D+01
  A( 7) =  -1.008066138D+01
  A( 8) =   6.010731909D-01
  A( 9) =  -1.567417031D-02 
  A=A[1:9]
  index=where(e ge 1.e3,nind)
  if nind gt 0 then begin
     sigma[1,index]=10.^(poly(alog10(e[index]),A))
  endif
  ;; fudge the cross sections for e < 1keV
  index=where(e lt 1.e3,nind)
  if nind gt 0 then begin
     sigma[1,index]=(e[index]*1.e-3)^.5*10.^(poly(alog10(1.e3),A))
  endif
  ;; m=3
  A=dblarr(20)
  A( 1) =   7.037287586D+04
  A( 2) =  -1.479161477D+05
  A( 3) =   1.370120708D+05 
  A( 4) =  -7.343180122D+04
  A( 5) =   2.509832081D+04
  A( 6) =  -5.674317075D+03
  A( 7) =   8.487767749D+02
  A( 8) =  -8.102284612D+01
  A( 9) =   4.480007503D+00 
  A(10) =  -1.093512342D-01
  A=A[1:10]
  index=where(e ge 2.e3,nind)
  if nind gt 0 then begin
     sigma[2,index]=10.^(poly(alog10(e[index]),A))
  endif
  ;; fudge the cross sections for e < 2keV
  index=where(e lt 2.e3,nind)
  if nind gt 0 then begin
     sigma[2,index]=(e[index]*1.e-3)^1.4*10.^(poly(alog10(2.e3),A))/2.8
  endif
  ;; m=4
  A=dblarr(20)
  A( 1) =   6.826447557D+04  
  A( 2) =  -1.431980004D+05
  A( 3) =   1.323968679D+05 
  A( 4) =  -7.083995050D+04
  A( 5) =   2.417608863D+04
  A( 6) =  -5.458418789D+03
  A( 7) =   8.154875237D+02 
  A( 8) =  -7.776012846D+01
  A( 9) =   4.295431731D+00
  A(10) =  -1.047567211D-01
  A=A[1:10]
  index=where(e ge 2.e3,nind)
  if nind gt 0 then begin
     sigma[3,index]=10.^(poly(alog10(e[index]),A))
  endif
  ;; fudge the cross sections for e < 2keV
  index=where(e lt 2.e3,nind)
  if nind gt 0 then begin
     sigma[3,index]=(e[index]*1.e-3)^2*10.^(poly(alog10(2.e3),A))/4.
  endif
  return,sigma
END
FUNCTION sigma2_adas,erel,mmax
  n=2.
  e=erel*1.e3  ;; conversion from [kev] to [eV]
  e=e*n^2      ;; correct for bug in the ADAS routines
  nel=n_elements(erel)
  sigma=replicate(0.,mmax,nel)
  index=where(e lt 1.e3,nind)
  if nind gt 0 then e[index]=1.e3
  ;; m=1
  ;sigma1=sigma1_adas(erel,mmax)
  sigma1=sigma1_adas(erel*n^2/1.^2,mmax)
  sigma[0,*]=sigma1[1,*]*1.^2/n^2;*n^4/1.^4
  ;; m=2
  As=dblarr(12)
  As( 1) =  -1.896015167D+06           
  As( 2) =   4.431727330D+06
  As( 3) =  -4.627815357D+06   
  As( 4) =   2.843068107D+06
  As( 5) =  -1.137952956D+06    
  As( 6) =   3.100801094D+05
  As( 7) =  -5.825744660D+04     
  As( 8) =   7.452319142D+03
  As( 9) =  -6.212350647D+02     
  As(10) =   3.047712749D+01
  As(11) =  -6.682658463D-01
  As=As[1:11] 
  sigma2s=10.^(poly(alog10(e),As))
  Ap=dblarr(12)
  Ap( 1) =  -1.614213508D+05          
  Ap( 2) =   3.772469288D+05
  Ap( 3) =  -3.924736424D+05  
  Ap( 4) =   2.393127027D+05
  Ap( 5) =  -9.470300966D+04  
  Ap( 6) =   2.541276100D+04
  Ap( 7) =  -4.682860453D+03  
  Ap( 8) =   5.851219013D+02
  Ap( 9) =  -4.744504549D+01  
  Ap(10) =   2.254460913D+00
  Ap(11) =  -4.767235839D-02
  Ap=Ap[1:11]
  sigma2p=10.^(poly(alog10(e),Ap))
  sigma[1,*]=0.25*sigma2s+0.75*sigma2p
  ;; fudge the m=2 cross section above 1.5e5
  index=where(erel gt 1.5e2,nind)
  if nind gt 0 then begin
     sigma2s=10.^(poly(alog10(1.5e5*4.),As))
     sigma2p=10.^(poly(alog10(1.5e5*4.),Ap))
     sigma[1,index]=(e[index]*1.e-3)^(-5.5)*(0.25*sigma2s+0.75*sigma2p)*2.e15
  endif
  ;;m=3
  As=dblarr(12)
  As( 1) =  -3.513030327D+05          
  As( 2) =   9.281116596D+05
  As( 3) =  -1.086843398D+06   
  As( 4) =   7.437325055D+05
  As( 5) =  -3.296609685D+05   
  As( 6) =   9.897503768D+04
  As( 7) =  -2.039707143D+04   
  As( 8) =   2.850670244D+03
  As( 9) =  -2.587092857D+02   
  As(10) =   1.377382945D+01
  As(11) =  -3.268306303D-01
  As=As[1:11] 
  sigma2s=10.^(poly(alog10(e),As))
  Ap=dblarr(12)
  Ap( 1) =  -1.901264631D+05         
  Ap( 2) =   5.124716103D+05
  Ap( 3) =  -6.101921504D+05    
  Ap( 4) =   4.234717934D+05
  Ap( 5) =  -1.899866398D+05    
  Ap( 6) =   5.764464326D+04
  Ap( 7) =  -1.199087959D+04     
  Ap( 8) =   1.689900512D+03
  Ap( 9) =  -1.545334374D+02     
  Ap(10) =   8.285001228D+00
  Ap(11) =  -1.978656474D-01
  Ap=Ap[1:11]
  sigma2p=10.^(poly(alog10(e),Ap))
  sigma[2,*]=0.25*sigma2s+0.75*sigma2p 

  ;; determine CX cross section in all m states from janev
  sigma_n=picx_janev(erel,n)-total(sigma,1) >0.
  ;; use scaling from LUIS 
  En=13.6/float(n)^2
  sigma_fudge=replicate(0.d0,mmax,nel)
  for m=0,mmax-1 do begin
     Em=13.6/float(m+1.)^2
     c=1.
     if total(sigma[m,*]) eq 0 then $
        sigma_fudge[m,*]=sigma_n/(c*sqrt(2.*!pi))*exp(-0.5*(En-Em)^2/c^2)
  endfor
  ;; normalize sigma to the JANEV cross sections
  sigma_tot=reform(total(sigma_fudge,1))
  factor=sigma_n/sigma_tot
  for m=0,mmax-1 do begin
     if total(sigma[m,*]) eq 0 then begin   
        sigma[m,*]=sigma_fudge[m,*]*factor
     endif
  endfor
  ;; if there are NANs in the array
  index=where(sigma ne sigma,nind)
  if nind gt 0 then sigma[index]=0.
  return,sigma
END




FUNCTION sigma3_adas,erel,mmax
  n=3.
  e=erel*1.e3  ;; conversion from [ev] to [keV]
  e=e*n^2      ;; correct for bug in the ADAS routines
  index=where(e lt 1.e3,nind)
  if nind gt 0 then e[index]=1.e3
  nel=n_elements(erel)
  sigma=replicate(0.,mmax,nel)
  ;; m=1
 ; sigma1=sigma1_adas(erel,mmax)
  sigma1=sigma1_adas(erel*n^2/1.^2,mmax)
  sigma[0,*]=sigma1[2,*]*1.^2/n^2;*n^4/1.^4
  ;;m=2
  A=dblarr(12)
  A( 1) =  -1.149224555D+06        
  A( 2) =   2.750368877D+06
  A( 3) =  -2.942222842D+06     
  A( 4) =   1.852584954D+06
  A( 5) =  -7.603284323D+05            
  A( 6) =   2.125284465D+05
  A( 7) =  -4.097580431D+04    
  A( 8) =   5.380901722D+03
  A( 9) =  -4.606297192D+02     
  A(10) =   2.321345254D+01
  A(11) =  -5.230186707D-01
  A=A[1:11] 
  sigma[1,*]=10^(poly(alog10(e),A)) 
 ;; fudge the m=2 cross section above 90keV
  index=where(erel gt 90.,nind)
  if nind gt 0 then begin
     sigma[1,index]=10.^(poly(alog10(9.e4*9.),A)) $
                    * (e[index]*1.e-3)^(-5.5)*1.e16
  endif 
  ;; m=3
  A=dblarr(11)
  A( 1) =  -4.302808608D+04          
  A( 2) =   9.499298161D+04
  A( 3) =  -9.264698488D+04     
  A( 4) =   5.236947172D+04
  A( 5) =  -1.890479538D+04     
  A( 6) =   4.519068626D+03
  A( 7) =  -7.152485009D+02     
  A( 8) =   7.227063167D+01
  A( 9) =  -4.230036444D+00    
  A(10) =   1.092702525D-01
  A=A[1:10] 
  sigma[2,*]=10^(poly(alog10(e),A)) 
  ;; fudge the m=3 cross section above 90keV
  index=where(erel gt 90.,nind)
  if nind gt 0 then begin
     sigma[2,index]=10.^(poly(alog10(9.e4*9.),A)) $
                    * (e[index]*1.e-3)^(-5.5)*.85e16
  endif

  ;; m=4
  A=dblarr(11)
  A( 1) =   1.705303425D+04  
  A( 2) =  -3.316878090D+04
  A( 3) =   2.792556433D+04    
  A( 4) =  -1.330264490D+04
  A( 5) =   3.921666688D+03   
  A( 6) =  -7.327555138D+02
  A( 7) =   8.476342861D+01  
  A( 8) =  -5.551987930D+00
  A( 9) =   1.577120745D-01
  A=A[1:9] 
  sigma[3,*]=10^(poly(alog10(e),A)) 
  ;; fudge the m=4 cross section above 90keV
  index=where(erel gt 90.,nind)
  if nind gt 0 then begin
     sigma[3,index]=10.^(poly(alog10(9.e4*9.),A)) $
                    * (e[index]*1.e-3)^(-5.5)*.82e16
  endif
  ;; m=5
  if mmax ge 5 then begin
     A=dblarr(12)
     A( 1) =  -2.786268232D+02        
     A( 2) =   4.269683825D+04
     A( 3) =  -8.973561028D+04  
     A( 4) =   8.365732310D+04
     A( 5) =  -4.524587937D+04  
     A( 6) =   1.563630402D+04
     A( 7) =  -3.580391824D+03 
     A( 8) =   5.432527332D+02
     A( 9) =  -5.267599631D+01  
     A(10) =   2.962329657D+00
     A(11) =  -7.362649692D-02
     A=A[1:11] 
     sigma[4,*]=10^(poly(alog10(e),A))
  endif
  
  if mmax ge 6 then begin
     ;;m=6
     ;; these are cross sections summed over all
     ;; higher states!@
     A=dblarr(12)
     A( 1) =   7.146969470D+05  
     A( 2) =  -1.665413326D+06
     A( 3) =   1.735840441D+06  
     A( 4) =  -1.065792786D+06
     A( 5) =   4.269334710D+05  
     A( 6) =  -1.165954977D+05
     A( 7) =   2.198700496D+04  
     A( 8) =  -2.827160468D+03
     A( 9) =   2.372409350D+02  
     A(10) =  -1.173264972D+01
     A(11) =   2.596865877D-01
     A=A[1:11]
     sigma_m6=10^(poly(alog10(e),A))  
     ;; fudge the m=6 cross section above 90keV
     index=where(erel gt 90.,nind)
     if nind gt 0 then begin
        sigma_m6[index]=10.^(poly(alog10(9.e4*9.),A)) $
                        * (e[index]*1.e-3)^(-7.)*2.e20
     endif 
     ;; use scaling from LUIS 
     En=13.6/float(n)^2
     sigma_fudge=replicate(0.d0,mmax,nel)
     for m=0,mmax-1 do begin
        Em=13.6/float(m+1.)^2
        c=1.
        if total(sigma[m,*]) eq 0 then $
           sigma_fudge[m,*]=sigma_m6/(c*sqrt(2.*!pi))*exp(-0.5*(En-Em)^2/c^2)
     endfor
     ;; normalize sigma to the JANEV cross sections
     sigma_tot=reform(total(sigma_fudge,1))
     factor=sigma_m6/sigma_tot
     for m=5,mmax-1 do begin
        if total(sigma[m,*]) eq 0 then begin  
           sigma[m,*]=sigma_fudge[m,*]*factor
        endif
     endfor
  endif
  ;; if there are NANs in the array
  index=where(sigma ne sigma,nind)
  if nind gt 0 then sigma[index]=0.
  return,sigma
END

FUNCTION sigman_fudge,erel,mmax,n
  nel=n_elements(erel)
  sigma=replicate(0.d0,mmax,nel)

  sigma1=sigma1_adas(erel*n^2/1.^2,mmax)
  sigma2=sigma2_adas(erel*n^2/2.^2,mmax)
  sigma3=sigma3_adas(erel*n^2/3.^2,mmax)
  ;; sigma1=sigma1_adas(erel,mmax)
  ;; sigma2=sigma2_adas(erel,mmax)
  ;; sigma3=sigma3_adas(erel,mmax)
  ;; for n=4:
  if n eq 4 then begin
     sigma[0,*]=sigma1[3,*]*1.^2/n^2;*n^4/1.^4
     sigma[1,*]=sigma2[3,*]*2.^2/n^2;*n^4/2.^4
     sigma[2,*]=sigma3[3,*]*3.^2/n^2;*n^4/3.^4
  endif
  ;; for n=5:
  if n eq 5 then begin
     sigma[1,*]=sigma2[4,*]*2.^2/n^2;*n^4/2.^4
     sigma[2,*]=sigma3[4,*]*3.^2/n^2;*n^4/3.^4
  endif
  ;; for n=6:
  if n eq 6 then begin
     sigma[1,*]=sigma2[5,*]*2.^2/n^2;*n^4/2.^4
     sigma[2,*]=sigma3[5,*]*3.^2/n^2;*n^4/3.^4
  endif  
  ;; determine CX cross section in all m states from janev
  sigma_n=picx_janev(erel,n)-total(sigma,1) >0.
  ;; use scaling from LUIS 
  En=13.6/float(n)^2
  sigma_fudge=replicate(0.d0,mmax,nel)
  for m=0,mmax-1 do begin
     Em=13.6/float(m+1.)^2
     c=1.
     if total(sigma[m,*]) eq 0 then $
        sigma_fudge[m,*]=sigma_n/(c*sqrt(2.*!pi))*exp(-0.5*(En-Em)^2/c^2)
  endfor
  ;; normalize sigma to the JANEV cross sections
  sigma_tot=reform(total(sigma_fudge,1))
  factor=sigma_n/sigma_tot
  for m=0,mmax-1 do begin
     if total(sigma[m,*]) eq 0 then begin   
        sigma[m,*]=sigma_fudge[m,*]*factor
     endif
  endfor
 ;; if there are NANs in the array
  index=where(sigma ne sigma,nind)
  if nind gt 1 then sigma[index]=0.
  if min(sigma) lt 0. then stop
  return,sigma
end
function sigma_cx,erel,n,m
;; charge exchange n and m resolved charge exchange cross-sections
;; H(+)+H(n)-->H(m)+H(+)
;; if m is set to -1 then return the janev Cross sections!
;; the n=4 state is derived by using the Principle of detailed balance
;; States higher than n=4 are taken from JANEV 
  if m eq -1 then return, picx_janev(erel,n)
  ;; cross sections from ADAS
  mmax=12
  case n of
     0: stop
     1:  sigma=sigma1_adas(erel,mmax)
     2:  sigma=sigma2_adas(erel,mmax)
     3:  sigma=sigma3_adas(erel,mmax)
     else: begin
        sigma=sigman_fudge(erel,mmax,n)
     end
  endcase
  sigma=reform(sigma[m-1,*])
  return, sigma
end


;; PROTON IMPACT IONIZATION
function omullane,energy,n
  ;; bgeiger proton impact ionization taken from the correction of
  ;; Omullane
  ;; the excited state ionization from janev seems to be wrong!
  ;; energy in [keV]
  n=float(n)
  if n le 5 then i=n-1. else i=4.
  Etil=energy*n^2  
;; n=     1,          2     ,     3     ,     4     ,     5
  b1= [2.0160e-03, 3.9330e-03, 1.1076e-02, 1.1033e-02, 1.1297e-02]
  b2= [3.7154e+00, 1.8188e+00, 1.6197e+00, 1.6281e+00, 1.8685e+00]
  b3= [3.9890e-02, 1.8870e-02, 6.7154e-03, 5.5955e-03, 1.5038e-02]
  b4= [3.1413e-01, 6.7489e-03, 5.1188e-03, 7.2023e-03, 1.1195e-01]
  b5= [2.1254e+00, 1.3768e+00, 1.8549e+00, 1.7358e+00, 1.0538e+00]
  b6= [6.3990e+03, 6.8852e+02, 2.3696e+02, 2.2755e+02, 8.6096e+02]
  b7= [6.1897e+01, 9.6435e+01, 7.8286e+01, 8.6339e+01, 8.9939e+01]
  b8= [9.2731e+03, 5.6515e+23, 1.0926e+23, 3.9151e+29, 1.9249e+04]    
  p1 = b1[i]*(i+1)^4
  p2 = Etil^b2[i] * exp(-b3[i]*Etil) / (1.d0 + b4[i]*Etil^b5[i])
  p3 = (b6[i]* exp(-b7[i]/Etil) *alog(1.d0  +b8[i]*Etil) ) /Etil
  sigma = 1.0d-16 * p1 * (p2 + p3)
  if n gt 5 then sigma=sigma*(n/5.)^4.
  return,sigma
END
;; PROTON IMPACT EXCTIATION
FUNCTION piexcitation_janev,eb,n,m
; D. Liu June 12,2007
; B. Geiger 2012
; H(+)+H(n)->H(+)+H(n'l')
; cross section for proton impact excitation
; in the mamal "Collision processes in low-temperature hydrogen plasmas"
; by Janev et al. 
; http://www.eirene.de/report_4105.pdf
; INPUT
;       eb      laboratory energy of proton (keV/amu)
;       n      lower energy state
;       m	upper energy state
; OUTPUT        sigma in cm^2
n=float(n)
m=float(m)
if n ge m then return,0.
if (n eq 1) then begin ;transition from the ground state
   if m eq 2 then begin
      a1=34.433
      a2=8.5476
      a3=7.8501
      a4=-9.2217
      a5=1.8020e-2
      a6=1.6931
      a7=1.9422e-3
      a8=2.9068
      a9=44.507
      a10=0.56870    
      sigma=a1*(a2*exp(-a3*eb)/eb^a4+a5*exp(-a6/eb)/(1.+a7*eb^a8)+$
                exp(-a9/eb)*alog(1.+a10*eb)/eb)*1.e-16
      return, sigma
   endif
   if m ge 3 then begin
      ;; n=1-> 3        4           5           6
      b1= [6.1950,     2.0661,     1.2449,     0.63771   ]
      b2= [5.5162e-3,  5.1335e-4,  3.0826e-4,  3.2949e-4 ]
      b3= [0.29114,    0.28953,    0.31063,    0.25757   ]
      b4= [-4.5264,   -2.2849,    -2.4161,    -2.2950    ]
      b5= [6.0311,     0.11528,    0.024664,   0.050796  ]
      b6= [-2.0679,   -4.8970,    -6.3726,    -5.5986    ]
      b7= [35.773,     34.975,     32.291,     37.174    ]
      b8= [0.54818,    0.91213,    0.21176,    0.39265   ]  
      if m le 6 then i=m-3 else i=3
      sigma=b1[i]*(b2[i]*exp(-b3[i]*eb)/$
                   (eb^b4[i]+b5[i]*eb^b6[i])+$
                   exp(-b7[i]/eb)*alog(1.+b8[i]*eb)/eb)*1.e-16  
      if m gt 6 then sigma = sigma * (6./m)^3
      return, sigma
   endif
endif
if n eq 2 then begin    	    
   ;; m    3       4         5
   c1=[394.51,   50.744,   18.264  ]
   c2=[0.013597, 0.014398, 0.013701]
   c3=[0.16565,  0.31584,  0.31711 ]
   c4=[-0.8949, -1.4799,  -1.4775  ]
   c5=[21.606,   19.416,   18.973  ]
   c6=[0.62426,  4.0262,   2.9056  ] 
   if m le 5 then i=m-3 else i=2
   sigma=c1[i]*(c2[i]*exp(-c3[i]*eb)/(eb^c4[i])+$
                exp(-c5[i]/eb)*alog(1.+c6[i]*eb)/eb)*1.e-16    
   if m gt 5 then begin
      A=[0.4610, 0.2475, 0.1465, 0.0920, 0.0605] 
      if m le 10 then ii=m-6 else ii=4 
      sigma=sigma*A[ii]
      if m gt 10 then sigma=sigma*(10./m)^3
   endif
   return, sigma
endif
if n eq 3 then begin
   ;; m 4         5          6
   c1=[1247.5,    190.59,   63.494    ]
   c2=[0.068781,  0.073307, 0.077953  ]
   c3=[0.521176,  0.54177,  0.53461   ]
   c4=[-1.2722,  -1.2894,  -1.2881    ]
   c5=[11.319,    11.096,   11.507    ]
   c6=[2.6235,    2.9098,   4.3417    ]
   if m le 6 then i=m-4 else i=2
   sigma=c1[i]*(c2[i]*exp(-c3[i]*eb)/(eb^c4[i])+$
                exp(-c5[i]/eb)*alog(1.+c6[i]*eb)/eb)*1.e-16
   if m gt 6 then begin
      A=[0.4670, 0.2545, 0.1540, 0.1000] 
      if m le 10 then ii=m-7 else ii=3
      sigma=sigma*A[ii]
      if m gt 10 then sigma=sigma*(10./m)^3
   endif
   return,sigma
endif  	      
if n ge 4 then begin
   n=float(n)
   m=float(m)
   etil=Eb/25.
   s=(m-n)
   D=exp(-1./(n*m*etil^2))
   A=8./(3.*s)*(m/(s*n))^3*(0.184-0.04/s^(2./3.))*(1.-0.2*s/(n*m))^(1.+2.*s)
   G=0.5*(etil*n^2./(m-1./m))^3.
   L=alog(1.+0.53*etil^2.*n*(m-2./m)/(1.+0.4*etil))
   F=(1.-0.3*s*D/(n*m))^(1.+2.*s)
     
   y=1./(1.-D*alog(18*s)/(4.*s))
   zpl=2./(etil*n^2*((2.-n^2/m^2)^0.5+1.))
   zmi=2./(etil*n^2*((2.-n^2/m^2)^0.5-1.))
   C2pl=zpl^2*alog(1.+2.*zpl/3.)/(2.*y+3.*zpl/2.)
   C2mi=zmi^2*alog(1.+2.*zmi/3.)/(2.*y+3.*zmi/2.)
   H=C2mi-C2pl
   sigma=8.8e-17*n^4/etil*(A*L*D+F*G*H)
   return, sigma
endif
end

;; --------------------------------
;; COLLISIONS WITH ELECTRONS:
;; --------------------------------
;; ELECTRON IMPACT IONIZATION
FUNCTION eiionization_janev2004,ecoll_keV,nint
  ;; INPUT
  ;;       ecoll      electron temperature (keV)
  ;;       n       energy level
  ;; OUTPUT        sigma in cm^2
  ecoll=ecoll_keV*1.e3
  n=float(nint)
  if n le 3 then begin
     An=replicate(0.d,6,3)
     ;; n=      1           2          3
     In=     [ 13.6,      13.6/2.^2, 13.6/3.^2]
     An[0,*]=[ 0.18450,   0.14784,   0.058463 ]
     An[1,*]=[-0.032226,  0.0080871,-0.051272 ]
     An[2,*]=[-0.034539, -0.062270,  0.85310  ]
     An[3,*]=[ 1.4003,    1.9414,   -0.57014  ] 
     An[4,*]=[-2.8115,   -2.1980,    0.76684  ]
     An[5,*]=[ 2.2986,    0.95894,   0.00     ]
     sum=0.d
     for j=1,5 do begin
        sum=sum+An[j,n-1]*(1.-In[n-1]/ecoll)^j  
     endfor
     sigma=1.e-13/(In[n-1]*ecoll)*(An[0,n-1]*alog(ecoll/In[n-1])+sum)   
     index=where(ecoll lt 2.,nind)
     if nind gt 0 then sigma[index]=0.
     return, sigma
  endif
  if n gt 3 then begin
     rn=1.94/n^1.57
     In=13.6/n^2
     xn=ecoll/In
     g0=0.9935+0.2328/n-0.1296/n^2
     g1=-1./n*(0.6282-0.5598/n+0.5299/n^2.)
     g2=1./n^2.*(0.3887-1.181/n+1.47/n^2.)
     An=32.*n/(3.*sqrt(3.)*!pi)*(g0/3.+g1/4.+g2/5.)
     b=1./n*(4.-18.63/n+36.24/n^2.-28.09/n^3.)
     Bn=2./3.*n^2.*(5.+b)
     sigma=1.76*n^2/xn*(1.-exp(-rn*xn))*$
           (An*alog(xn)+(Bn-An*alog(2.*n^2))*(1.-1./xn)^2)*1.e-16
     index=where(xn le 1.,nind)
     if nind gt 0 then sigma[index] = 0.
     return, sigma
  endif
end
;; ELECTRON IMPACT EXCITATION
function fnm, nint,mint
  n=float(nint)
  m=float(mint)
  x=1.-(n/m)^2.
  if n eq 1 then g=[1.133,-0.4059,0.0714]
  if n eq 2 then g=[1.0785,-0.2319,0.02947]
  if n ge 3 then begin
     g=fltarr(3)
     g[0]=0.9935 + 0.2328/n - 0.1296/n^2
     g[1]=-1./n * (0.6282 - 0.5598/n + 0.5299/n^2)
     g[2]= 1./n^2*(0.3887 - 1.1810/n + 1.4700/n^2)
  endif   
  g=g[0]+g[1]/x+g[2]/x^2
  fnm=32./(3.*sqrt(3.)*!pi) * n/m^3* 1/x^3 *g
  return, fnm
end
FUNCTION eiexcitation_janev2004,ecoll_keV,nint,mint
  ;; INPUT
  ;;       Energy      electron temperature (keV)
  ;;       n      lower energy state
  ;;       m	upper energy state
  ;; OUTPUT        sigma in cm^2
  ecoll=ecoll_keV*1.e3
  n=float(nint)
  m=float(mint)
  sigma=fltarr(n_elements(ecoll))
  if n ge m then stop
  if (n eq 1) then begin        ;transition from the ground state
     if m eq 2 then begin
        sigma0=5.984
        deltaE=10.2
        a=0.228 
        b=0.1865 
        c=0.5025
        x=ecoll/deltaE
        An=[4.4979,1.4182,-20.877,49.735,-46.249,17.442]
        index=where(ecoll gt 10.2 and ecoll le 11.56,nind)
        if nind gt 0 then sigma[index]=(a+b*(ecoll[index]-deltaE))*1.e-16
        index=where(ecoll ge 11.56 and ecoll le 12.23,nind)
        if nind gt 0 then sigma[index] = replicate(0.,nind)+c*1.e-16
        index=where(ecoll le 10.2 or ecoll gt 12.23,nind)
        if nind gt 0 then begin
           sum=0.0d
           for j=1,5 do begin
              sum=sum+An[j]/x[index]^(j-1.)
           endfor	
           sigma[index]=sigma0/(deltaE*x[index])*(An[0]*alog(x[index])+sum) $
                        *1.e-16
        endif
        index=where(x le 1.,nind)
        if nind gt 0 then sigma[index] = 0.
        return, sigma
     endif
     if m ge 3. and m le 5. then begin
        sigma0=5.984
        ;; m     3         4        5   
        deltaE= [12.09,   12.75,   13.06]
        alpha=  [0.38277, 0.41844, 0.45929]
        An=replicate(0.d,5,3)
        An[0,*]=[0.75448, 0.24300, 0.11508]
        An[1,*]=[0.42956, 0.24846, 0.13092]
        An[2,*]=[-0.58288,0.19701, 0.23581]
        An[3,*]=[1.0693,  0.00,    0.00]
        An[4,*]=[0.00,    0.00,    0.00]
        x=ecoll/deltaE[m-3]
        sum=0.d
        for j=1,4 do begin
           sum=sum+An[j,m-3]/x^(j-1.)
        endfor 
        sigma=sigma0/(deltaE[m-3]*x)*(1.-1./x)^alpha[m-3]*$
              (An[0,m-3]*alog(x)+sum)*1.e-16  
        index=where(x le 1.,nind)
        if nind gt 0 then sigma[index] = 0.
        return, sigma 
     endif
     if m gt 5. then begin
        deltaEm=13.6*(1.-1./m^2)
        xm=ecoll/deltaEm
        ym=1.-(1./m)^2
        r=0.45
        f1m=fnm(1.,m)
        Am=2.*n^2*f1m/ym
        Bm=4./(m^3*ym)*(1.+4./(3.*ym)-0.603/ym^2)
        sigma=1.76e-16/(ym*xm)*(1.-exp(-r*ym*xm))*$
              (Am*(alog(xm)+1./(2.*xm))+(Bm-Am*alog(2./ym)) $
               *(1.-1./xm))
        index=where(xm le 1.,nind)
        if nind gt 0 then sigma[index] = 0.
        return,sigma
     endif
  endif
  ;; transition from excited states
  if n eq 2 and m eq 3 then begin
     sigma0=5.984
     alpha=1.3196 
     An=[38.906, 5.2373, 119.25, -595.39, 816.71]
     deltaE=13.6*(1./n^2-1./m^2.)
     x=ecoll/deltaE
     sum=0.d
     for j=1,4 do begin
        sum=sum+An[j]/x^(j-1)
     endfor      	 
     sigma=sigma0/(deltaE*x)*(1.-1./x)^alpha*(An[0]*alog(x)+sum)*1.e-16
     index=where(x le 1.,nind)
     if nind gt 0 then sigma[index] = 0.
     return,sigma
  endif
  ;; all the other transitions
  deltaEnm=13.6*(1./n^2-1./m^2)
  xnm=ecoll/deltaEnm
  ynm=1.-(n/m)^2
  rn=1.94/n^1.57
  fnm=fnm(n,m)
  Anm=2.*n^2*fnm/ynm
  bn=1./n*(4.-18.63/n+36.24/n^2-28.09/n^3)
  Bnm=4.*n^4/(m^3*ynm^2)*(1.+4./(3.*ynm)+ bn/ynm^2.)     
  sigma=1.76e-16*n^2/(ynm*xnm)*(1.-exp(-rn*ynm*xnm))*$
        (Anm*(alog(xnm)+1./(2.*xnm))+(Bnm-Anm*alog(2.*n^2./ynm)) $
         *(1.-1./xnm))
  index=where(xnm le 1.,nind)
  if nind gt 0 then sigma[index] = 0.
  return,sigma
end
  





;; --------------------------------
;; COLLISIONS WITH IMPURITIES
;; --------------------------------
;; CHARGE EXCHANGE
FUNCTION icx_adas,erel,nint,q
  sigma=replicate(0.,n_elements(erel))
  e=erel*1.e3
  ;; below 10 eV, set energy to 10eV!
  index=where(e lt 10.,nind)
  if nind gt 0 then e[index]=10.
  A=replicate(0.d0,12)
  n=float(nint)
  if q eq 5 then begin
     case n of
        1:begin
           index=where(e lt 1.e3,nind)
           if nind gt 0 then e[index]=1.e3
           A( 1) =   1.174052518D+03         
           A( 2) =  -1.793561728D+03
           A( 3) =   1.117522436D+03      
           A( 4) =  -3.679435571D+02
           A( 5) =   6.750816878D+01    
           A( 6) =  -6.542029074D+00
           A( 7) =   2.614113716D-01
           A=A[1:10]
           sigma=10^(poly(alog10(e),A))
           index=where(e gt 400.e3,nind)
           if nind gt 0 then sigma[index]=0.
           return, sigma
        end
        2:begin
           A( 1) =   6.603246818D+01    
           A( 2) =  -3.072575676D+02
           A( 3) =   5.030801019D+02    
           A( 4) =  -4.585636345D+02
           A( 5) =   2.568666393D+02    
           A( 6) =  -9.185150382D+01
           A( 7) =   2.100012584D+01   
           A( 8) =  -2.964174788D+00
           A( 9) =   2.346396110D-01       
           A(10) =  -7.943766873D-03
           A=A[1:10]
           sigma=10^(poly(alog10(e),A))
           index=where(erel gt 1.e2,nind)
           if nind gt 0 then sigma[index]=0.
           return,sigma
        end
        else: return,sigma
     endcase
  endif
  if q eq 6 then begin
     index=where(e lt 1.5e3,nind)
     if nind gt 0 then e[index]=1.5e3
     case n of
        1:begin
           A( 1) =   2.007882674D+02      
           A( 2) =  -3.546893286D+02
           A( 3) =   2.381542403D+02     
           A( 4) =  -8.355431742D+01
           A( 5) =   1.617519888D+01    
           A( 6) =  -1.638152470D+00
           A( 7) =   6.768953863D-02
        end
        2:begin
           A( 1) =   9.151879441D+05        
           A( 2) =  -2.134573133D+06
           A( 3) =   2.223792624D+06    
           A( 4) =  -1.362648703D+06
           A( 5) =   5.438401343D+05   
           A( 6) =  -1.477110500D+05
           A( 7) =   2.764972254D+04    
           A( 8) =  -3.522105245D+03
           A( 9) =   2.921934171D+02    
           A(10) =  -1.425552507D+01
           A(11) =   3.106007048D-01
        end
        3: begin
           A( 1) =   9.208877916D+05  
           A( 2) =  -2.147294379D+06
           A( 3) =   2.236451628D+06    
           A( 4) =  -1.370042347D+06
           A( 5) =   5.466461899D+05      
           A( 6) =  -1.484338816D+05
           A( 7) =   2.777765778D+04      
           A( 8) =  -3.537459450D+03
           A( 9) =   2.933884362D+02      
           A(10) =  -1.430994136D+01
           A(11) =   3.117002878D-01
        end
        else: return,sigma
     endcase
     A=A[1:11]
     e=e*n^2
     sigma=10^(poly(alog10(e),A))
  endif
  return,sigma
END
FUNCTION icx_janev,erel,nint,qint
  ;; INPUTS
  ;; e : energy [keV/amu]
  ;; n : state 
  ;; q : impurity charge
  ;; OUTPUT
  ;; sigma [cm^2]
  e=erel
  n=float(nint)
  q=float(qint)
 ;; IF BORON then specific data
  if n eq 1 and q eq 5 then begin
  ;; page 166 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
  ;; Volume 4 (by Janev and Smith 1993
     A=fltarr(12)
     A[1]  = 31.226
     A[2]  = 1.1442
     A[3]  = 4.8372e-8
     A[4]  = 3.0961e-10
     A[5]  = 4.7205
     A[6]  = 6.2844e-7
     A[7]  = 3.1297
     A[8]  = 0.12556
     A[9]  = 0.30098
     A[10] = 5.9607e-2
     A[11] =-0.57923
     sigma=1.e-16*A[1]*(exp(-A[2]/e^A[8]) $
                        /(1.+A[3]*e^2+A[4]*e^A[5]+A[6]*e^A[7]) $
                        +A[9]*exp(-A[10]*e)/e^A[11])
  endif 
 ;; IF CARBON then specific data
  if n eq 1 and q eq 6 then begin
  ;; page 168 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
  ;; Volume 4 (by Janev and Smith 1993
     A=fltarr(12)
     A[1]  = 418.18
     A[2]  = 2.1585
     A[3]  = 3.4808e-4
     A[4]  = 5.3333e-9
     A[5]  =4.6556
     A[6]  =0.33755
     A[7]  = 0.81736
     A[8]  = 0.27874
     A[9]  = 1.8003e-6
     A[10] = 7.1033e-2
     A[11] =0.53261
     sigma=1.e-16*A[1]*(exp(-A[2]/e^A[8]) $
                        /(1.+A[3]*e^2+A[4]*e^A[5]+A[6]*e^A[7]) $
                        +A[9]*exp(-A[10]*e)/e^A[11])
  endif   
  ;; General data
 if n gt 1 or (q ne 5 and q ne 6) then begin
     ;; page 174 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
     ;; Volume 4 (by Janev and Smith 1993
    A=1.507e5
    B=1.974e-5
    etilde=e*n^2/q^0.5
    sigma=q*n^4*7.04e-16*A/(etilde^3.5*(1+B*etilde^2))* $
          (1.-exp(-2.*etilde^3.5*(1.+B*etilde^2)/(3*A)))
  endif
  return,sigma
end
FUNCTION icx_adas_janev,e,n,q
  ;; INPUTS
  ;; e : energy [keV/amu]
  ;; n : state 
  ;; q : impurity charge
  ;; OUTPUT
  ;; sigma [cm^2]
  sigma = icx_adas(e,n,q)
  if total(sigma) le 0. then sigma = icx_janev(e,n,q)
  return,sigma
end
;; IMPURITY IMPACT IONIZATION
FUNCTION iiionization,erel,nint,qint
  ;; INPUTS
  ;; e : energy [keV/amu]
  ;; n : state 
  ;; q : impurity charge
  ;; OUTPUT
  ;; sigma [cm^2]
  e=erel
  n=float(nint)
  q=float(qint)
  ;; IF BORON then specific data
  if n eq 1 and q eq 5 then begin
;; page 152 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
;; Volume 4 (by Janev and Smith 1993
     A=fltarr(9)
     A[1] = 351.52
     A[2] = 233.63
     A[3] = 3.2952e3
     A[4] = 5.3787e-6
     A[5] = 1.8834E-2
     A[6] =-2.2064
     A[7] = 7.2074
     A[8] =-3.78664
     sigma=1.e-16*A[1]*(exp(-A[2]/e)*alog(1+A[3]*e)/E $
                        +A[4]*exp(-A[5]*e)/((e^A[6])+A[7]*(e^A[8])))
  endif
  ;; IF CARBON then specific data!
  if n eq 1 and q eq 6 then begin
;; page 154 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
;; Volume 4 (by Janev and Smith 1993
     A=fltarr(9)
     A[1] = 438.36
     A[2] = 327.10
     A[3] = 1.4444e5
     A[4] = 3.5212e-3
     A[5] = 8.3031e-3
     A[6] =-0.63731
     A[7] = 1.9116e4
     A[8] =-3.1003
     sigma=1.e-16*A[1]*(exp(-A[2]/e)*alog(1+A[3]*e)/E $
                        +A[4]*exp(-A[5]*e)/((e^A[6])+A[7]*(e^A[8])))
  endif
  ;; General data
  if n gt 1 or (q ne 5 and q ne 6) then begin
     ;; page 160 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
     ;; Volume 4 (by Janev and Smith 1993
     M=sqrt(0.283)
     B=4.04
     gamma=0.662
     c=137.
     v=sqrt(e/25.)
     u=n*v
     sigma_b=3.52e-16*n^4*q^2/u^2*(M^2*(alog(u^2/(c^2-u^2))-u^2/c^2) $
                                   +B-gamma/u^2) > 0.
     lambda=0.76d0
     sigma=exp(-lambda*q/u^2)*sigma_b
  endif
  return,sigma
end
;; IMPURITY IMPACT EXCITATION
FUNCTION iiexcitation,erel,nint,mint,qint
  ;; INPUTS
  ;; e : energy [keV/amu]
  ;; n : state 
  ;; q : impurity charge
  ;; OUTPUT
  ;; sigma [cm^2]
  n=float(nint)
  m=float(mint)
  q=float(qint)
  e=erel
;;General data
;;n=1:
  if n eq 1 then begin
     ;; page 132-134 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
     ;; Volume 4 (by Janev and Smith 1993
     if m eq 2 then begin
        A=fltarr(7)
        A[1]  = 38.738
        A[2]  = 37.033
        A[3]  = 0.39862
        A[4]  = 7.7582e-5
        A[5]  = 0.25402
        A[6]  =-2.7418
     endif
     if m eq 3 then begin
        A=fltarr(7)
        A[1]  = 4.3619
        A[2]  = 57.451
        A[3]  = 21.001
        A[4]  = 2.3292e-4
        A[5]  = 0.083130
        A[6]  =-2.2364
     endif
     if m eq 4 then begin
        A=fltarr(7)
        A[1]  = 1.3730
        A[2]  = 60.710
        A[3]  = 31.797
        A[4]  = 2.0207e-4
        A[5]  = 0.082513
        A[6]  =-2.3055
     endif
     if m ge 5 then begin
        A=fltarr(7)
        A[1]  = 0.56565
        A[2]  = 67.333
        A[3]  = 55.290
        A[4]  = 2.1595e-4
        A[5]  = 0.081624
        A[6]  = -2.1971
     endif  
     etilde=e/q
     xsi=2.^(0.5238*(1-sqrt(2./q)))
     sigma=q*1.e-16*xsi*A[1]*(exp(-A[2]/etilde)*alog(1+A[3]*etilde)/etilde $
                              + A[4]*exp(-A[5]*etilde)/etilde^A[6])
     if m gt 5 then sigma=sigma*(5./m)^3
     return,sigma
  endif
;;n=2:
  if n eq 2 then begin
     ;; page 138 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
     ;; Volume 4 (by Janev and Smith 1993
     if m eq 3 then begin
        A=fltarr(7)
        A[1]  = 358.03
        A[2]  = 25.283
        A[3]  = 1.4726
        A[4]  = 0.014398
        A[5]  = 0.12207
        A[6]  =-0.86210
     endif
     if m eq 4 then begin
        A=fltarr(7)
        A[1]  = 50.744
        A[2]  = 19.416
        A[3]  = 4.0262
        A[4]  = 0.014398
        A[5]  = 0.31584
        A[6]  =-1.4799
     endif
     if m ge 5 then begin
        A=fltarr(7)
        A[1]  = 18.264
        A[2]  = 18.973
        A[3]  = 2.9056
        A[4]  = 0.013701
        A[5]  = 0.31711
        A[6]  = -1.4775
     endif
     etilde=e/q
     xsi=2.^(0.5238*(1-sqrt(2./q)))
     sigma=q*1.e-16*xsi*A[1]*(exp(-A[2]/etilde)*alog(1+A[3]*etilde)/etilde $
                              + A[4]*exp(-A[5]*etilde)/etilde^A[6])
     if m gt 5 then begin
        hi=2.^(0.397*(1.-sqrt(2./q)))
        A=[0.4610, 0.2475,0.1465,0.092, 0.0605]
        if m le 10 then ii=m-6 else ii=4
        sigma=sigma*hi*A[ii]
        if m gt 10 then sigma=sigma*(10./m)^3
     endif
     return,sigma
  endif
;;n=3:
  if n eq 3 then begin
     ;; page 142 in ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FUSION
     ;; Volume 4 (by Janev and Smith 1993
     if m eq 4 then begin
        A=fltarr(7)
        A[1]  = 1247.5
        A[2]  = 11.319
        A[3]  = 2.6235
        A[4]  =0.068781
        A[5]  = 0.521176
        A[6]  =-1.2722
     endif
     if m eq 5 then begin
        A=fltarr(7)
        A[1]  = 190.59
        A[2]  = 11.096
        A[3]  = 2.9098
        A[4]  = 0.073307
        A[5]  = 0.54177
        A[6]  =-1.2894
     endif
     if m ge 6 then begin
        A=fltarr(7)
        A[1]  = 63.494
        A[2]  = 11.507
        A[3]  = 4.3417
        A[4]  = 0.077953
        A[5]  = 0.53461
        A[6]  =-1.2881
     endif
     etilde=e/q
     xsi=2.^(0.397*(1-sqrt(2./q)))
     sigma=q*1.e-16*xsi*A[1]*(exp(-A[2]/etilde)*alog(1+A[3]*etilde)/etilde $
                              + A[4]*exp(-A[5]*etilde)/etilde^A[6])
     if m gt 6 then begin
        hi=2.^(0.397*(1.-sqrt(2./q)))
        A=[0.4670, 0.2545,0.1540,0.1]
        if m le 10 then ii=m-7 else ii=3
        sigma=sigma*hi*A[ii]
        if m gt 10 then sigma=sigma*(10./m)^3
     endif
     return,sigma
  endif
  if n ge 4 then begin
     n=float(n)
     m=float(m)
     etil=E/(25.*q)
     hi=2.^(0.322*(1.-sqrt(2./q)))
     s=(m-n)
     D=exp(-1./(n*m*etil^2))
     A=8./(3.*s)*(m/(s*n))^3*(0.184-0.04/s^(2./3.))*(1.-0.2*s/(n*m))^(1.+2.*s)
     G=0.5*(etil*n^2./(m-1./m))^3.
     L=alog(1.+0.53*etil^2.*n*(m-2./m)/(1.+0.4*etil))
     F=(1.-0.3*s*D/(n*m))^(1.+2.*s)
     
     y=1./(1.-D*alog(18*s)/(4.*s))
     zpl=2./(etil*n^2*(sqrt(2.-n^2/m^2)+1.))
     zmi=2./(etil*n^2*(sqrt(2.-n^2/m^2)-1.))
     C2pl=zpl^2*alog(1.+2.*zpl/3.)/(2.*y+3.*zpl/2.)
     C2mi=zmi^2*alog(1.+2.*zmi/3.)/(2.*y+3.*zmi/2.)
     H=C2mi-C2pl
     sigma=q*hi*8.86e-17*n^4/etil*(A*D*L+F*G*H)
     return, sigma
  endif
end
