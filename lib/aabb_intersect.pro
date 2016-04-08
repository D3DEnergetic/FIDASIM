PRO aabb_intersect, rc, dr, r0, d0, intersect, r_enter, r_exit
    ;+##`aabb_intersect, rc, dr, r0, d0, length, ri, rf`
    ;+Calculates intersection length of a ray and an axis aligned bounding box (AABB)
    ;+###Input Arguments
    ;+     **rc**: Center of AABB
    ;+
    ;+     **dr**: [length, width, height] of AABB
    ;+
    ;+     **r0**: starting point of ray
    ;+
    ;+     **d0**: direction of ray
    ;+
    ;+###Output Arguments
    ;+     **intersect**: Intersection length of ray and AABB
    ;+
    ;+     **ri**: Optional, ray enterence point
    ;+
    ;+     **rf**: Optional, ray exit point
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> aabb_intersect, [0,0,0], [1,1,1], [-1,0,0], [1,0,0], intersect, ri, rf
    ;+IDL> print, intersect
    ;+    1.0
    ;+IDL> print, ri
    ;+    -0.5  0.0  0.0
    ;+IDL> print, rf
    ;+     0.5  0.0  0.0
    ;+```

    v0 = d0/sqrt(total(d0*d0))

    ;; There are 6 sides to a cube/grid
    side_inter = dblarr(6)
    ;; Intersection points of ray with planes defined by grid
    ipnts = dblarr(3,6)

    ;; Find whether ray intersects each side
    for i=0L,5 do begin
        j = fix(floor(i/2))
        ind = where([0,1,2] ne j)
        if abs(v0[j]) gt 0 then begin
            ;; Intersection point with plane
            ipnts[*,i] = r0 + v0*( ( (rc[j] + ( (i mod 2)-0.5)*dr[j] ) - r0[j])/v0[j] )
            ;; Check if point on plane is within grid side
            if abs(ipnts[ind[0],i] - rc[ind[0]]) le 0.5*dr[ind[0]] and $
               abs(ipnts[ind[1],i] - rc[ind[1]]) le 0.5*dr[ind[1]] then side_inter[i]=1
        endif
    endfor

    intersect = 0.0
    r_enter = r0
    r_exit = r0
    w = where(side_inter ne 0,nw)
    if nw ge 2 then begin
        ;;Find two unique intersection points
        nunique = 0
        for i=0,nw-2 do begin
            if total(ipnts[*,w[0]] eq ipnts[*,w[i+1]]) ne 3 then begin
                w = [w[0],w[i+1]]
                nunique = 2
                break
            end
        end

        if nunique eq 2 then begin
            vi = ipnts[*,w[1]]-ipnts[*,w[0]]
            vi = vi/sqrt(total(vi*vi))
            dot_prod = total(v0*vi)
            if dot_prod gt 0.0 then begin
                r_enter = ipnts[*,w[0]]
                r_exit = ipnts[*,w[1]]
            endif else begin
                r_enter = ipnts[*,w[1]]
                r_exit = ipnts[*,w[0]]
            endelse
            ;; Calculate intersection length
            intersect = sqrt(total((r_exit-r_enter)^2.0))
        endif
    endif
END

