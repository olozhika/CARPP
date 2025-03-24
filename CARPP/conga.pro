    function conga,a,fwhm

; Convolve a square image "a" with Gaussian of specified fwhm [pixels].  

; Embed into an even-dimensioned array if necessary.
    nx = n_elements(a(*,0))
    ny = n_elements(a(0,*))
    if nx ne ny then begin
	message,/cont,'array must be square'
	return,-1
    endif

    if nx mod 2 eq 1 then begin
	b = fltarr(nx+1, nx+1) ;even
	b(1:nx,1:nx) = a
    endif else b = a  ;even dimention

; Transform to spatial frequency domain.
    sff = fft(b)

; Apply Gaussian taper.
    N = n_elements(b(*,0))
    sff = shift(sff, N/2, N/2)
    w = 4.*alog(2.)*N/(!pi*fwhm)
    r = shift(dist(N,N),N/2,N/2)
    arg = -4.*alog(2.)*(r/w)^2
    gauss = fltarr(N,N)
    ok = where(arg gt -20.)
    gauss(ok) = exp(arg(ok))
    sff = sff*gauss
    sff = shift(sff, N/2, N/2)

; Transform back again.
    b = float(fft(sff, /inverse))
    if nx mod 2 eq 1 then b = b(1:nx, 1:nx)
    return,b

    end
