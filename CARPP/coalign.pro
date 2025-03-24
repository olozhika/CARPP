    function coalign,psf,refpsf,xyshift,flux,FITWIDTH=fitwidth,NARRATE=narrate

; Shift a PSF based on a least squares fit with a reference PSF.
;
; K. A. Marsh 11/25/03; modified 7/21/06
;
; Input parameters:
;	psf	= input PSF array
;	refpsf	= reference PSF array
;
; Output parameters:
;	xyshift	= 2-element array containing the applied x and y shifts [pixels]
;	flux	= estimated flux, such that <coaligned psf> = flux*refpsf + nu
;
; Optional keywords:
;	FITWIDTH= width of fitting region [pixels]. Default: use entire image.
;	NARRATE	= if set, type out the applied (x,y) shift [pixels]
;
; Returned function value: Shifted (coaligned) PSF array.

    common sumsqcom,a,r,trs,as,alpha

    ftol = 1.e-4		; tolerance for minimizations

    if keyword_set(fitwidth) then begin
	ix = n_elements(psf(*,0))/2
	iy = n_elements(psf(0,*))/2
	h = (fitwidth/2) < (ix < iy)
	a = psf(ix-h:ix+h-1, iy-h:iy+h-1)
	r = refpsf(ix-h:ix+h-1, iy-h:iy+h-1)
    endif else begin
	a = psf
	r = refpsf
    endelse

    trs = total(r^2)
    xyshift = [0.,0.]
    xi = transpose([[1.0, 0.0],[0.0, 1.0]])
    powell,xyshift,xi,ftol,fmin,'sumsq'
    if keyword_set(narrate) then print,'    Shift x,y: ',xyshift
    flux = alpha

    if keyword_set(fitwidth) then begin
	x0 = xyshift(0)
	y0 = xyshift(1)
	as = shiftrot(psf,x0,y0)
    endif

    return,as

    end
