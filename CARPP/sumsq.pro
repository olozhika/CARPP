    function sumsq,p

; Calculate the sum of squares of residuals for the maximum likelihood solution
; of a shifted-psf fit.  (Called by COALIGN).
;
; K. A. Marsh, 11/25/03

    common sumsqcom

    x0 = p(0)
    y0 = p(1)
    as = shiftrot(a,x0,y0)
    alpha = total(as*r)/trs
    return, -total(as*r)

    end
