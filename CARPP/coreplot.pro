    pro coreplot,infiles,outname,SNR=snr,RCHISQMAX=rchisqmax

; Plot results of COREFIT.

    if not keyword_set(rchisqmax) then rchisqmax = 1.e35
    nfiles = n_elements(infiles)
    mass = fltarr(nfiles)
    beta1 = fltarr(nfiles)

    for n = 0,nfiles-1 do begin
	h = headfits(infiles(n))
	mass(n) = sxpar(h,'COREMASS')
	beta1(n) = sxpar(h,'BETA')
    endfor

    if nfiles ge 4 then begin
	x = moment(mass)
	cmass = x(0)
	sigcmass = sqrt(x(1))
	x = moment(beta1)
	cbeta1 = x(0)
	sigcbeta1 = sqrt(x(1))
    endif else begin
	cmass = total(mass)/nfiles
	cbeta1 = total(beta1)/nfiles
	sigcmass = 0.
	sigcbeta1 = 0.
    endelse

    profile0 = readfits(infiles(0),h,/silent)
    npoints = sxpar(h,'NAXIS1')
    cell = sxpar(h,'CDELT2')*3600.

    !p.multi = [0,1,2]
    device = 'ps'
    set_plot, device
    outfile = outname + '.ps'
    print,'Writing out file '+outfile
    device, file=outfile

    r = findgen(npoints)*cell
    title = outname
    yrange = [0.,1.5]
    if keyword_set(snr) then title = title + $
	string(format='("   SNR =",f6.1)',snr)
    plot,r,profile0(*,1),yrange=yrange,/ystyle,xtitle='r [arcsec]', $
	ytitle='Relative density',title=title
    if nfiles gt 1 then begin
	for n = 1,nfiles-1 do begin
	    profile = readfits(infiles(n),h,/silent)
	    if sxpar(h,'RCHISQ') le rchisqmax then oplot,r,profile(*,1),line=1
	endfor
	xyouts,0.45*(npoints-1)*cell,0.9*max(yrange), $
	  string(format='("m_cloud(est) =",f6.2," +/-",f6.2," msun")',$
	  cmass,sigcmass)
	xyouts,0.45*(npoints-1)*cell,0.8*max(yrange), $
	  string(format='("beta(est)    =",f6.2," +/-",f6.2)',$
	  cbeta1,sigcbeta1)
    endif else begin
	xyouts,0.45*(npoints-1)*cell,0.9*max(yrange), $
	  string(format='("m_cloud(est) =",f6.2," msun")',cmass)
	xyouts,0.45*(npoints-1)*cell,0.8*max(yrange), $
	  string(format='("beta(est)    =",f6.2)',cbeta1)
    endelse

    ymax = -1.e35
    for n = 0,nfiles-1 do begin
	profile = readfits(infiles(n),h,/silent)
	if max(profile(*,0)) gt ymax then ymax = max(profile(*,0))
    endfor
    yrange = fltarr(2)
    yrange(1) = 1.1*ymax
    plot,r,profile0(*,0),yrange=yrange,/ystyle,xtitle='r [arcsec]', $
	ytitle='T [K]'
    if nfiles gt 1 then begin
	for n = 1,nfiles-1 do begin
	    profile = readfits(infiles(n),h,/silent)
	    if sxpar(h,'RCHISQ') le rchisqmax then oplot,r,profile(*,0),line=1
	endfor
    endif
    device, /close

    return

    end

