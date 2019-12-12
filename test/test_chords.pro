FUNCTION test_chords
    
    ;; Chords
    ulens = [170.d0,170.d0,170.d0]
    vlens = dblarr(3)
    wlens = replicate(150.d0,3)
    lens = transpose([[ulens],[vlens],[wlens]])
    ulos = [200.d0,170.d0,140.d0]
    vlos = dblarr(3)
    wlos = dblarr(3)
    los = transpose([[ulos],[vlos],[wlos]])
    axis = los - lens
    for i=0,2 do begin
        axis[*,i] = axis[*,i]/sqrt(total(axis[*,i]^2.0))
    endfor
    
    sigma_pi = replicate(1.d0,3)
    spot_size = replicate(0.d0,3)
    radius = sqrt(ulos^2.d0 + vlos^2.d0) 
    id = ["f1","f2","f3"]
    chords = {nchan:3L,system:"SPECTRAL",id:id,data_source:'test_chords.pro',$
              lens:lens,axis:axis,spot_size:spot_size, sigma_pi:sigma_pi, $
              radius:radius}
    return, chords
END
