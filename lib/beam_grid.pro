FUNCTION beam_grid, nbi, rstart, $
    nx = nx, ny = ny, nz = nz, dv=dv, $
    length=length, width=width, height=height

    beam_grid = {err:1}

    if not keyword_set(length) then length = 100.0 ;cm
    if not keyword_set(width) then width = 80.0 ;cm
    width = max(width,2*nbi.widy)
    if not keyword_set(height) then height = 80.0 ;cm
    width = max(height,2*nbi.widz)

    if not keyword_set(dv) then dv = 8.0 ;cm^3

    dv3 = dv^(1.0/3.0)
    if not keyword_set(nx) then nx = round(length/dv3)
    if not keyword_set(ny) then ny = round(width/dv3)
    if not keyword_set(nz) then nz = round(height/dv3)

    xmin = 0.d0
    xmax = double(length)
    ymin = -width/2.d0
    ymax = width/2.d0
    zmin = -height/2.d0
    zmax = height/2.d0

    src = nbi.src
    axis = nbi.axis/sqrt(total(nbi.axis^2))
    pos = src + 100*axis
    rsrc = sqrt(src[0]^2 + src[1]^2)
    if sqrt(src[0]^2 + src[1]^2) lt rstart then begin
        error, "Source radius cannot be less then rstart"
        goto, GET_OUT
    endif

    dis = sqrt(total((src - pos)^2.0))
    beta = double(asin((src[2]-pos[2])/dis))
    alpha = double(atan((pos[1]-src[1]),(pos[0]-src[0])))
    gamma = 0.d0 
    a = axis[0]^2 + axis[1]^2
    b = 2*(src[0]*axis[0] + src[1]*axis[1])
    c = src[0]^2 + src[1]^2 - rstart^2
    t = (-b - sqrt(b^2 - 4*a*c))/(2*a)
    origin = src + t*axis

    beam_grid = {nx:nx, ny:ny, nz:nz, $
                 xmin:xmin, xmax:xmax, $
                 ymin:ymin, ymax:ymax, $
                 zmin:zmin, zmax:zmax, $
                 alpha:alpha,beta:beta, gamma:gamma, $
                 origin:origin }

    GET_OUT:
    return, beam_grid
    
END
