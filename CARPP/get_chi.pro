function get_chi, imagefiles,imagesigma,wavelen,fwhm,back,ncells,grainsize,    $
	distance,outname,SEDFILE=sedfile,RZERO=rzero,INDSUM=indsum,           $
	RCHISQMAX=rchisqmax,MODELIMAGES=modelimages,APERTUREMAP=aperturemap,  $
	SEDPLOT=sedplot,XPEAKPIX=xpeakpix,YPEAKPIX=ypeakpix,ADDNOISE=addnoise,$
	ROUTSET=routset,ALPHA_PRE=alpha_pre,BETA_PRE=beta_pre,INPUT_P=input_p,$
	FIT_SETTING=fit_setting,START_PIX=start_pix,quanzhong=quanzhong

    common corecom,nlayers,rout,a,d,imageset,nxset,nyset,lambda1set,           $
	fwhmset,varset,npixels,phimin,iterations,finished,spectra,sigspectra, $
	lambda1spec,isum,sedmask,coverage,modsed,obsapflux,pwrlaw,rz,use_sed, $
	gammainit,alpha_preset,beta_preset
    common dencom,npeak,pc_for_pixel,phicore_quanzhong,full_sphere
    common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, $
	arcsec, cloud_d
    common p_for_phicore,p_fixed,p_cloudmass

    a = grainsize
    d = distance
    setup_header,a=a,d=d
    fit_rout=0 ;bad result, don't use
    
    ; renew central pix coordinate from possible DS9 method to
    ; normal array that starts from [0,0]
    if keyword_set(xpeakpix) then begin
      if not keyword_set(start_pix) then begin
        start_pix=1
        ;print,'Asm peakpixs start from [1,1] as in DS9, Plz set start_pix=0 if wrong'
      endif
      xpeakpix=xpeakpix-start_pix
      ypeakpix=ypeakpix-start_pix
    endif
    
    ; set fitting methods
    if not keyword_set(fit_setting) then fit_setting=[0,0,0,0,1]
    if fit_rout then fit_setting=[0,0,0,0,0]
    if keyword_set(quanzhong) then phicore_quanzhong=quanzhong

    gammainit=2.;original, it is 2. Tpf index
    nlayers = 50		; assumed number of layers in model
    ;if nlayers gt 50 then nlayers = 50
    
    if keyword_set(alpha_pre) then alpha_preset=alpha_pre
    if keyword_set(beta_pre) then beta_preset=beta_pre
    
    if not keyword_set(rchisqmax) then rchisqmax = 1.5
    if not keyword_set(spectitle) then spectitle = 'Observed SED of core'
    pwrlaw = keyword_set(rzero)
    if pwrlaw then rz = rzero
    if keyword_set(indsum) then isum = indsum else isum = -1.e35

; Preprocess input data. This includes regridding the images.
    prefit,imagefiles,imagesigma,wavelen,fwhm,back,ncells,racent,deccent, $
	error,SEDFILE=sedfile,XPEAKPIX=xpeakpix,YPEAKPIX=ypeakpix,        $
	ADDNOISE=addnoise
    if error then return,-1
    infiles=imagefiles
    for temp=0,n_elements(infiles)-1 do begin
      infiles[temp]='regrid_'+infiles[temp]
    endfor
    ;infiles = ['regrid_','regrid_','regrid_']+imagefiles
    use_sed = keyword_set(sedfile)
    if use_sed then inspec = sedfile
    outlabel = outname

    if keyword_set(addnoise) then begin
	infiles = 'noisy_'+infiles
	if use_sed then inspec = 'noisy_'+sedfile
	label = string(format='("sample",i3.3,"_")',addnoise)
	outlabel = label+outname
    endif

; Read the regridded images back in.
    nbands = n_elements(infiles)
    lambda1set = fltarr(nbands)
    nxset = intarr(nbands)
    nyset = intarr(nbands)
    fwhmset = fltarr(nbands)
    varset = fltarr(nbands)
    nxmax = 0
    nymax = 0
    for n = 0,nbands-1 do begin
	head = headfits(infiles(n))
	nxset(n) = sxpar(head,'NAXIS1')
	nyset(n) = sxpar(head,'NAXIS2')
	if nxset(n) gt nxmax then nxmax = nxset(n)
	if nyset(n) gt nymax then nymax = nyset(n)
	lambda1set(n) = sxpar(head,'WAVELEN')
	;fwhmset(n) = sxpar(head,'FWHM');in arcsec
	fwhmset(n) = sxpar(head,'FWHM')/(sxpar(head,'CDELT2')*3600.);in pixel, changed by XYC,bcz in conga what we need is in pix
	varset(n) = sxpar(head,'SIGOBS')^2 ;sigobs is assumed uncertainy depending on the input imagesigma
        ;SEEHERE ^2
    endfor
    pixel = sxpar(head,'CDELT2')
    npixels = sxpar(head,'NPIXELS')
    if keyword_set(ROUTSET) then rout=ROUTSET else $ ;added by XYC
      rout = 0.5*npixels*pixel*!dtor*d	; assume that the
					; original image size corresponds to
					; the outer edge of the cloud
					; rout in [pc]
		if keyword_set(ROUTSET) then full_sphere=0 else full_sphere=1
					
    npixels = round(npixels)
    pc_for_pixel=pixel*!dtor*d

    imageset = fltarr(nxmax,nymax,nbands)
    for n = 0,nbands-1 do imageset(0:nxset(n)-1, 0:nyset(n)-1, n) = $
;	readfits(infiles(n),/silent) ;commented out on Oct 25, 2012 due to
;/silent syntax error
;  try without option

	readfits(infiles(n))

; Read in spectra.
    if use_sed then sedread, inspec,aperture,slitwid,pa_slit,ra_ap,dec_ap, $
	lambda1spec,spectra,sigspectra,Dspec

; Create a set of masks representing the responses of SED apertures.
  if use_sed then begin
    nap = n_elements(ra_ap)
    nlambda1 = n_elements(lambda1spec)
    sedmask = fltarr(npixels,npixels,nap,nlambda1)
    nbig = 9*npixels
    off = 4*npixels
    apmap = fltarr(nbig,nbig)
    for iap = 0,nap-1 do begin
	mask = fltarr(nbig,nbig)
	crval1 = sxpar(head,'CRVAL1')
	crval2 = sxpar(head,'CRVAL2')
	cd = cos(crval2*!dtor)
	x0 = sxpar(head,'CRPIX1') - 1 - (ra_ap(iap)-crval1)*cd/pixel
	y0 = sxpar(head,'CRPIX2') - 1 + (dec_ap(iap)-crval2)/pixel
	ilo = (round(x0+off - 0.5*slitwid/(pixel*3600.))>0)<(nbig-1)
	ihi = round(x0+off + 0.5*slitwid/(pixel*3600.)) < (nbig-1)
	jlo = (round(y0+off - 0.5*aperture/(pixel*3600.))>0)<(nbig-1)
	jhi = round(y0+off + 0.5*aperture/(pixel*3600.)) < (nbig-1)
	mask(ilo:ihi, jlo:jhi) = 1.
	mask = rot(mask, pa_slit, 1., x0+off, y0+off, /pivot, $
	    cubic=-0.5)
	apmap = apmap + mask
	apfunc = mask(off:off+npixels-1,off:off+npixels-1)
	for n = 0,nlambda1-1 do begin
	    psfwid = ((lambda1spec(n)*1.e-6/Dspec)/!dtor)/pixel
	    apresp = conga(apfunc,psfwid)
	    sedmask(*,*,iap,n) = conga(apresp,psfwid)
	endfor
    endfor
    if keyword_set(aperturemap) then begin
	apmap(nbig/2-3:nbig/2+3,nbig/2) = 2.
	apmap(nbig/2,nbig/2-3:nbig/2+3) = 2.
	apmap(nbig/2-npixels/2:nbig/2+npixels/2, nbig/2-npixels/2) = 2.
	apmap(nbig/2-npixels/2:nbig/2+npixels/2, nbig/2+npixels/2) = 2.
	apmap(nbig/2-npixels/2, nbig/2-npixels/2:nbig/2+npixels/2) = 2.
	apmap(nbig/2+npixels/2, nbig/2-npixels/2:nbig/2+npixels/2) = 2.
	apmapfile = outlabel+'_apertures.fits'
	;print,'Writing out '+apmapfile
	writefits,apmapfile,apmap
    endif
    modsed = fltarr(nlambda1,nap)
    obsapflux = fltarr(nbands,nap)
  endif ; if use_sed

; Run a grid of models to get initial parameters for Newton-Raphson step.
; Assume the core is optically thin at the longest wavelength.
    longest = max(lambda1set,nlw)
    ok = where(imageset eq imageset and imageset gt 0, count); original its imageset ne 0, changed by XYC
    if count eq 0 then begin
       message,/CONT,'No valid data'
       return,-1
    endif
    gseed = long(1.e9*(imageset(ok(0))- min(imageset(ok)))/ $
            max(imageset(ok)) - min(imageset(ok)))
    
    if 1 then begin
      imageset=imageset>0.
      ;index_array=fltarr(npixels,npixels)
      cent_pix=(npixels-1)/2
      for indx_i = 0, npixels-1 do begin
        for indx_j = 0, npixels-1 do begin
          dist1 = sqrt((indx_i-cent_pix)^2+(indx_j-cent_pix)^2) > 1.
          if (dist1 gt cent_pix) then imageset[indx_i,indx_j,*]=0
        endfor
      endfor
    endif
   
    maxpasses = 5
    notconverged = 1
    ipass = 0
    
    nstate = 7
    if pwrlaw then nstate = nstate -1		; number of unknowns
    if keyword_set(alpha_preset) then nstate = nstate -1
    if keyword_set(beta_preset) then nstate = nstate -1
    
    if keyword_set(alpha_preset) then alpha_i=1 else alpha_i=0
    if keyword_set(beta_preset) then beta_i=1 else beta_i=0
    
    phimin = 1.e35
    chimin = 1.e35
    iterations = 0
    finished = 0

    if keyword_set(input_p) then p=input_p ;XYC's input p
    
    if fit_setting(0) gt 0 then chi = phicore_noqz(p[1:*])
    if fit_setting(4) gt 0 then chi = phicore(p)
    return,chi

    end
