    pro ploterrorbars,x,y,sigx,sigy,MAGNITUDES=magnitudes

; Plot error bars in both axes.  If the inputs represent magnitudes, then the
; MAGNITUDES keyword should be set.

    N = n_elements(x)
    xx = fltarr(2)
    yy = fltarr(2)

    if keyword_set(magnitudes) then begin
	tle = 2.5*alog10(exp(1.))
	for i = 0,N-1 do begin
	    xx(0) = x(i) - 2.5*alog10(1. + sigx(i)/tle)
	    xx(1) = x(i) - 2.5*alog10((1. - sigx(i)/tle)>1.e-35)
	    yy(*) = y(i)
	    oplot,xx,yy
	    yy(0) = y(i) - 2.5*alog10(1. + sigy(i)/tle)
	    yy(1) = y(i) - 2.5*alog10((1. - sigy(i)/tle)>1.e-35)
	    xx(*) = x(i)
	    oplot,xx,yy
	endfor
    endif else begin
	for i = 0,N-1 do begin
	    xx(0) = x(i) - sigx(i)
	    xx(1) = x(i) + sigx(i)
	    yy(*) = y(i)
	    oplot,xx,yy
	    yy(0) = y(i) - sigy(i)
	    yy(1) = y(i) + sigy(i)
	    xx(*) = x(i)
	    oplot,xx,yy
	endfor
    endelse

    return

    end
