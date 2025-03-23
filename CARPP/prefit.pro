    pro prefit,imagefiles,imagesigma,wavelen,fwhm_beam,fwhm_resolu,back,ncells,racent,deccent,$
	error,SEDFILE=sedfile,XPEAKPIX=xpeakpix,YPEAKPIX=ypeakpix,           $
	ADDNOISE=addnoise, throw_value=throw_value

; Prepare images and SEDs for processing with COREFIT.
;
; Input parameters:
;	imagefiles	=	name of FITS files of input images
;	imagesigma	=	fractional uncertainty in observed images, with
;				respect to peak
;	wavelen		= 	wavelengths of input images [microns]
;	fwhm		=	FWHMs of PSFs [arcsec]
;	back		=	background levels [Jy/beam]
;	ncells		=	size of model image in units of pixel size of
;				first image
; Output parameters:
;	racent		=	RA of center of model image [deg] 
;	deccent		=	Dec of center of model image [deg] 
;
; Optional keywords:
;	SEDFILE		=	name of text file containing SEDs
;	ADDNOISE	=	if set (to an integer value between 1 and 999),
;				add simulated Gaussian noise to measurements
;	XPEAKPIX	=	x coordinates [pixels] of core centers
;	YPEAKPIX	=	y coordinates [pixels] of core centers

    error = 0
    nimages = n_elements(imagefiles)
    if ncells mod 2 eq 0 then ncells = ncells - 1 ; ne 0, ncells +1, changed by XYC
    icent = (ncells-1)/2 ; ncells/2; changed by XYC
    h = headfits(imagefiles(0))
    cdelt1 = sxpar(h,'CDELT1')
    cdelt2 = sxpar(h,'CDELT2')
    crval1 = sxpar(h,'CRVAL1')
    crval2 = sxpar(h,'CRVAL2')
    crpix1 = sxpar(h,'CRPIX1')
    crpix2 = sxpar(h,'CRPIX2')
    if cdelt1 eq 0 or cdelt2 eq 0 or crval1 eq 0 or crval2 eq 0 or $
	crpix1 eq 0 or crpix2 eq 0 then begin
	print,'One or more of CDELT, CRVAL, CRPIX keywords missing from ' $
	    + imagefiles(0)
	error = 1
	return
    endif
    pixel = cdelt2*3600.
    r = shift(dist(ncells,ncells),icent,icent)
    out = where(r gt icent)
    beam = fwhm_resolu(0)/pixel
    ref = exp(-4.*alog(2.)*(r/beam)^2)
    if keyword_set(addnoise) then seed = addnoise
    print,'Reading input images ...'

    for n = 0,nimages-1 do begin
;	a = readfits(imagefiles(n),h,/silent)
	a = readfits(imagefiles(n),h)
	if sxpar(h,'CDELT2') eq 0 then begin
	    print,'No CDELT2 present in '+imagefiles(n)
	    error = 1
	    return
	endif
	nx = sxpar(h,'NAXIS1')
	;peak = max(a,loc)
	;ipeak = loc mod nx
	;jpeak = loc/nx
	if keyword_set(xpeakpix) and keyword_set(ypeakpix) then begin
	    ipeak = round(xpeakpix(n));ipeak = round(xpeakpix(n) - 1.) the -1 is in other place
	    jpeak = round(ypeakpix(n));jpeak = round(ypeakpix(n) - 1.)
	    peak = a(ipeak,jpeak)
	    register = 0
	endif else begin
	    peak = max(a,loc)
	    ipeak = loc mod nx
	    jpeak = loc/nx
	    register = 1
	endelse
	print,imagefiles(n),':   peak at (',ipeak,',',jpeak,')'
	
	; make sure that the image covers the ncells ;UNFINISHED, by XYC
	;size_temp=size(a,/Dimensions)
	;if_ncells_ok=(max(where(ipeak-icent lt 0)) eq -1)* $
  ;             (max(where(ipeak+icent gt size_temp[0]-1)) eq -1)* $
	;             (max(where(jpeak-icent lt 0)) eq -1)* $
  ;             (max(where(jpeak+icent gt size_temp[1]-1)) eq -1)
	;if if_ncells_ok ne 1 then
	
	mag = sxpar(h,'CDELT2')*3600./pixel ;20230227 suofang beishu
	
	;extend a  added by XYC
	a_todo=[ipeak-icent,sxpar(h,'NAXIS1')-(ipeak+icent),$
	  jpeak-icent,sxpar(h,'NAXIS2')-(jpeak+icent)]
	a_if= a_todo ge 0
	b=a((ipeak-icent)>0:(ipeak+icent)<(sxpar(h,'NAXIS1')-1),$
	  (jpeak-icent)>0:(jpeak+icent)<(sxpar(h,'NAXIS2')-1)) - back(n)
	
	if total(a_if) lt 4 then begin
	  b1=b
	  b=fltarr(ncells,ncells)
	  b(abs(0<a_todo[0]):ncells-1-abs(0<a_todo[1]),$
	    abs(0<a_todo[2]):ncells-1-abs(0<a_todo[3]))=b1
	endif
		
	;b = a(ipeak-icent:ipeak+icent,jpeak-icent:jpeak+icent) - back(n)
  ;changed by XYC, it was b = a(ipeak-icent:ipeak+icent-1, jpeak-icent:jpeak+icent-1) - back(n)
	if n ne 0 then begin
	    b = rot(b,0.,mag,icent,icent,/pivot,cubic=-0.5);suofang, mag is the ratio
	    sxaddpar,h,'CDELT1',-pixel/3600.
	    sxaddpar,h,'CDELT2',pixel/3600.
	endif
	if register then begin
	    b = coalign(b,ref,xyshift,flux)
	    print,'    Applied shift [pixels]:',xyshift
	endif else xyshift = fltarr(2)

	if n eq 0 then begin
	    cd = cos(crval2*!dtor)
	    racent = crval1 + (ipeak-crpix1 - xyshift(0))*cdelt1/cd ;Changed by XYC, it was racent = crval1 + (ipeak-crpix1+1 - xyshift(0))*cdelt1/cd
	    deccent = crval2 + (jpeak-crpix2 - xyshift(1))*cdelt2; same as above
	endif

  if keyword_set(throw_value) then b[where(b lt throw_value)] = 0

; Convert intensity units from Jy/beam to Jy/pixel
    beam = fwhm_beam(n)/pixel
    ref = exp(-4.*alog(2.)*(r/beam)^2) ; XYC doesnt understand what is r doing here, but the ans is correct
    ref = ref/total(ref)
    barea = 1./max(ref)
    b = b/barea

	sxaddpar,h,'CDELT1',cdelt1
	sxaddpar,h,'CDELT2',cdelt2
	sxaddpar,h,'CRPIX1',icent
	sxaddpar,h,'CRPIX2',icent
	sxaddpar,h,'CRVAL1',racent
	sxaddpar,h,'CRVAL2',deccent
	sxaddpar,h,'BUNIT','Jy/pixel' ;added by XYC
	sxaddpar,h,'WAVELEN',wavelen(n)
	sxaddpar,h,'FWHM',fwhm_resolu(n)
        if n_elements(imagesigma) eq 1 then $
	  sigobs = imagesigma*max(b) else sigobs = imagesigma(n)*max(b)
        
	sxaddpar,h,'SIGOBS',sigobs ;SIGOBS is a set-ted value, not a robust value
	sxaddpar,h,'NPIXELS',ncells;ncells-1; XYC doesnt understand why there is a minus1, so she removed it
	image = b ;image = b(1:ncells-1,1:ncells-1); XYC doesnt understand why throw the edges, so she removed it
	outimage = 'regrid_'+imagefiles(n)
	if keyword_set(addnoise) then begin
	    seed = seed + 1000
	    image = image + sigobs*randomn(seed,ncells,ncells);randomn(seed,ncells-1,ncells-1) ;XYC doesnt understand why there is a minus1, so she removed it
	    outimage = 'noisy_' + outimage
	endif
	print,'    Writing out regridded image: '+outimage
	print,'    Total flux =',total(image),' Jy'
	writefits,outimage,image,h
    endfor

    if keyword_set(addnoise) and keyword_set(sedfile) then begin
	sedread, sedfile,aperture,slitwid,pa_slit,ra_ap,dec_ap,lambda1spec, $
	    spectra,sigspectra
	nlambda1 = n_elements(lambda1spec)
	nap = n_elements(ra_ap)
	seed = seed + 1000
	spectra = spectra + sigspectra*randomn(seed,nlambda1,nap)
	openr,1,sedfile
	openw,2,'noisy_'+sedfile
	print,'Writing out noisy SED: noisy_'+sedfile
	line = ' '
	done = 0
	while not done do begin
	    readf,1,line
	    printf,2,line
	    if strpos(line,'Jy') ne -1 then done = 1
	endwhile
	close,1
	nfc = 2*nap
	outform = string(format='("(f6.3,",i2,"f11.5)")',nfc)
	for i = 0,nlambda1-1 do begin
	    fluxcols = fltarr(nfc)
	    fluxcols(2*indgen(nap)) = spectra(i,*)
	    fluxcols(2*indgen(nap)+1) = sigspectra(i,*)
	    printf,2,format=outform,lambda1spec(i),fluxcols
	endfor
	close,2
    endif

    return

    end
