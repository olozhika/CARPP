    function shiftrot,array,xmove,ymove
;
; Shift an image by (xmove,ymove) pixels, using IDL ROT function.
;
    nx = n_elements(array(*,0))
    ny = n_elements(array(0,*))
    if round(abs(xmove)) ge nx/2 or round(abs(ymove)) ge ny/2 then $
	return, 0.*array

    xcent = (nx - 1.)/2.
    ycent = (ny - 1.)/2.
    arrayshift = rot(array,0.,1.,xcent-xmove,ycent-ymove,cubic=-0.5)

; Suppress wraparound effect.
    big = fltarr(2*nx,2*ny)
    big(nx-nx/2:nx+nx/2, ny-ny/2:ny+ny/2) = arrayshift ;fixed by XYC, it was big(nx-nx/2:nx+nx/2-1, ny-ny/2:ny+ny/2-1) = arrayshift
    xnew = nx + round(xmove)
    ynew = ny + round(ymove)
    mask = big
    mask(xnew-nx/2:xnew+nx/2, ynew-ny/2:ynew+ny/2) = !values.f_nan
    out = where(mask eq mask)
    big(out) = 0
    arrayshift = big(nx-nx/2:nx+nx/2, ny-ny/2:ny+ny/2);fixed by XYC, it was xxx-1
    return, arrayshift

    end
