FUNCTION beam_grid, nbi, rstart, $
    nx = nx, ny = ny, nz = nz, dv=dv, $
    length=length, width=width, height=height
    ;+#beam_grid
    ;+ Calculates settings for a grid that aligns with the neutral beam.
    ;+***
    ;+##Arguments
    ;+    **nbi**: Neutral beam geometry structure
    ;+
    ;+    **rstart**: Radial start position of beam grid [cm]
    ;+
    ;+##Keyword Arguments
    ;+    **dV**: Cell volume [\(cm^3\)]: Defaults to 8.0
    ;+
    ;+    **nx**: Number of cells in length: Default determined by `dV`
    ;+
    ;+    **ny**: Number of cells in width: Default determined by `dV`
    ;+
    ;+    **nz**: Number of cells in height: Default determined by `dV`
    ;+
    ;+    **length**: Length of grid along beam sightline. [cm]: Defaults to 100 cm
    ;+
    ;+    **width**: Width of grid [cm]: Defaults to 50 cm
    ;+
    ;+    **height**: Height of grid [cm]: Defaults to 50 cm
    ;+
    ;+##Return Value
    ;+    Structure containing beam grid settings suitable for the [Namelist File](|url|/page/03_technical/01_input_files/01_namelist_file.html)
    ;+   
    ;+##Example Usage
    ;+```idl
    ;+IDL> grid = beam_grid(nbi,200.0,nx=100,ny=50,nz=50,length=100,width=50,height=50)
    ;+```


    beam_grid = {err:1}

    if not keyword_set(length) then length = 100.0 ;cm

    if not keyword_set(width) then width = 80.0 ;cm
    if width lt nbi.widy then  begin
        warning, "Grid width is smaller then the source width"
        print, "width: ", width
        print, "source width: ",nbi.widy
    endif
 
    if not keyword_set(height) then height = 80.0 ;cm
    if height lt nbi.widz then begin
        warning, "Grid height is smaller then the source height" 
        print, "height: ", height
        print, "source height: ",nbi.widz
    endif

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

    beam_grid = {nx:fix(nx), ny:fix(ny), nz:fix(nz), $
                 xmin:xmin, xmax:xmax, $
                 ymin:ymin, ymax:ymax, $
                 zmin:zmin, zmax:zmax, $
                 alpha:alpha,beta:beta, gamma:gamma, $
                 origin:origin }

    GET_OUT:
    return, beam_grid
    
END
