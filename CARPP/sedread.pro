    pro sedread, sedfile,aperture,slitwid,pa_slit,ra_ap,dec_ap,lambda1spec, $
	spectra,sigspectra,Dspec

; Read in SED data.

    Dspec = 0.85	; telescope aperture for SED observations [m]
    slitwid = 19.6	; slit width [arcsec]
    apra = $
      ['ra_1', 'ra_2', 'ra_3', 'ra_4', 'ra_5', 'ra_6', 'ra_7', 'ra_8', 'ra_9']
    apdec = $
      ['dec_1','dec_2','dec_3','dec_4','dec_5','dec_6','dec_7','dec_8','dec_9']
    naplist = n_elements(apra)
    ra_ap = dblarr(naplist)
    dec_ap = dblarr(naplist)
    lambda1spec = fltarr(10000)
    lambda1 = 0.

    line = ' '
    openr,1,sedfile
    nap = 0
    i = -1

    while not eof(1) do begin
	readf,1,line
	parts = str_sep(line,'#')
	if n_elements(parts) gt 1 and strpos(line,'source') eq -1 $
	  and strpos(line,'=') ne -1 then begin
	    mark = strpos(parts(1),'=')
	    value = 0.d0
	    reads,strmid(parts(1),mark+1,100),value
	    if strpos(line,'# aperture') ne -1 then aperture = value
	    if strpos(line,'# cdelt') ne -1 then cdelt = value
	    if strpos(line,'# pa_slit') ne -1 then pa_slit = value
	    for n = 0,naplist-1 do begin
		if strpos(line,'# '+apra(n)) ne -1 then begin
		    nap = nap + 1
		    ra_ap(nap-1) = value
		endif
		if strpos(line,'# '+apdec(n)) ne -1 then $
		    dec_ap(nap-1) = value
	    endfor
	endif

	if strpos(line,'Jy') ne -1 then begin
	    fluxcols = fltarr(2*nap)
	    spectra = fltarr(10000,nap)
	    sigspectra = fltarr(10000,nap)
	endif

	if strpos(line,'#') eq -1 then begin
	    reads,line,lambda1,fluxcols
	    i = i+1
	    lambda1spec(i) = lambda1
	    spectra(i,*) = fluxcols(2*indgen(nap))
	    sigspectra(i,*) = fluxcols(2*indgen(nap)+1)
	endif
    endwhile
    close,1
    aperture = aperture*cdelt
    ra_ap = ra_ap(0:nap-1)
    dec_ap = dec_ap(0:nap-1)
    lambda1spec = lambda1spec(0:i)
    spectra = spectra(0:i,*)
    sigspectra = sigspectra(0:i,*)
    return

    end
