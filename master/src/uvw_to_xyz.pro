FUNCTION uvw_to_xyz, alpha, beta, gamma, uvw, origin=origin
    ;+##`uvw_to_xyz(alpha, beta, gamma, uvw, origin=[0,0,0])`
    ;+ Express non-rotated coordinate `uvw` in rotated `xyz` coordinates
    ;+###Arguments
    ;+     **alpha**: Rotation angle about z [radians]
    ;+
    ;+     **beta**: Rotation angle about y' [radians]
    ;+
    ;+     **gamma**: Rotation angle about x" [radians]
    ;+
    ;+     **xyz**: Point in rotated coordinate system
    ;+
    ;+###Keyword Arguments
    ;+     **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates.
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> xyz = uvw_to_xyz(!DPI/2,0.0,!DPI/3,uvw)
    ;+```
    if not keyword_set(origin) then origin=[0.0,0.0,0.0]
    s = size(uvw,/dim)
    if n_elements(s) ne 2 then s=[s,1]
    uvw_shifted = transpose(uvw - tile_array(origin,1,s[1]))

    R = transpose(tb_zyx(alpha,beta,gamma))

    xyz = R##uvw_shifted
    
    return, transpose(xyz)
END
