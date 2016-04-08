FUNCTION xyz_to_uvw, alpha, beta, gamma, xyz, origin = origin
    ;+##`xyz_to_uvw(alpha, beta, gamma, xyz, origin=[0,0,0])`
    ;+Express rotated coordinate `xyz` in non-rotated `uvw` coordinates
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
    ;+IDL> uvw = xyz_to_uvw(!DPI/2,0.0,!DPI/3,xyz)
    ;+```
    if not keyword_set(origin) then origin=[0.0,0.0,0.0]
    s = size(xyz,/dim)
    if n_elements(s) ne 2 then s=[s,1]
    xyz = transpose(xyz) ;Column Vector

    R = tb_zyx(alpha,beta,gamma)

    uvw = R##xyz 

    return, tranpose(uvw) + tile_array(origin,1,s[1])
END

