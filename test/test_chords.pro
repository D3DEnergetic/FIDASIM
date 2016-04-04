FUNCTION test_chords
    
    ;; Chords
    ulens = dblarr(3)
    vlens = [-170.d0,-170.d0,-170.d0]
    wlens = replicate(100.d0,3)
    lens = transpose([[ulens],[vlens],[wlens]])
    ulos = dblarr(3)
    vlos = [-200.d0,-170.d0,-140.d0]
    wlos = dblarr(3)
    los = transpose([[ulos],[vlos],[wlos]])
    axis = los - lens
    for i=0,2 do begin
        axis[*,i] = axis[*,i]/sqrt(total(axis[*,i]^2.0))
    endfor
    
    sigma_pi = replicate(1.d0,3)
    spot_size = replicate(0.d0,3)
    radius = sqrt(ulos^2.d0 + vlos^2.d0) 
    chords = {nchan:3,system:"SPECTRAL",data_source:'test_chords.pro',$
              lens:lens,axis:axis,spot_size:spot_size, sigma_pi:sigma_pi, $
              radius:radius}
    return, chords
END
