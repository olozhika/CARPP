    pro carpp, imagefiles,imagesigma,wavelen,fwhm_beam,back,ncells,grainsize,$
	distance,outname,SEDFILE=sedfile,RZERO=rzero,INDSUM=indsum,           $
	RCHISQMAX=rchisqmax,MODELIMAGES=modelimages,APERTUREMAP=aperturemap,  $
	SEDPLOT=sedplot,XPEAKPIX_IN=xpeakpix_in,YPEAKPIX_IN=ypeakpix_in,ADDNOISE=addnoise,$
	ROUTSET=routset,ALPHA_PRE=alpha_pre,BETA_PRE=beta_pre,INPUT_P=input_p,$
	FIT_SETTING=fit_setting,START_PIX=start_pix,weight_index=weight_index,      $
  gamma_pre=gamma_pre,double=double,fwhm_resolu=fwhm_resolu,fit_rout=fit_rout, $
  makeround=makeround, precise_pre=precise_pre,throw_value=throw_value, $
  MAX_RESOLVE_SIZE=max_resolve_size

    common corecom,nlayers,rout,a,d,imageset,nxset,nyset,lambda1set,           $
	fwhmset,varset,npixels,phimin,iterations,finished,spectra,sigspectra, $
	lambda1spec,isum,sedmask,coverage,modsed,obsapflux,pwrlaw,rz,use_sed, $
	gammainit,alpha_preset,beta_preset,max_resolve,printf_p_init
    common dencom,npeak,pc_for_pixel,phicore_quanzhong,full_sphere
    common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, $
	arcsec, cloud_d
    common p_for_phicore,p_fixed,p_cloudmass
    common parameter_range, delta_beta, delta_alpha, param_Tmin, delta_cloudmass

    a = grainsize
    d = distance
    setup_header,a=a,d=d
    setup_parameter_range
    
    if not keyword_set(MAX_RESOLVE_SIZE) then max_resolve=fltarr(n_elements(imagefiles)) $
      else max_resolve=MAX_RESOLVE_SIZE
    ;else max_resolve will be like [0.,0.,10.] in arcsec 

    if not keyword_set(precise_pre) then precise=1 else precise=precise_pre
    if not keyword_set(fwhm_resolu) then fwhm_resolu=fwhm_beam

    ; renew central pix coordinate from possible DS9 method to
    ; normal array that starts from [0,0]
    if keyword_set(xpeakpix_in) then begin
      if not keyword_set(start_pix) then begin
        start_pix=1
        print,'Asm peakpixs start from [1,1] as in DS9, Plz set start_pix=0 if wrong'
      endif
      xpeakpix=xpeakpix_in-start_pix
      ypeakpix=ypeakpix_in-start_pix
    endif
    
    ; set fitting methods
    if not keyword_set(fit_setting) then fit_setting=[0,-1]
    if keyword_set(fit_rout) then fit_setting=[0,0]
    if keyword_set(weight_index) then phicore_quanzhong=weight_index

    print,'-----------------------------------------'
    print,'Begin COREFIT on '+systime()
    print,'-----------------------------------------'
    
    nlayers = 50		; assumed number of layers in model
    ;if nlayers lt (ncells-1)/2 then nlayers = (ncells-1)/2
    
    if keyword_set(alpha_pre) then alpha_preset=alpha_pre
    if keyword_set(beta_pre) then beta_preset=beta_pre
    if keyword_set(gamma_pre) then gammainit=gamma_pre else gammainit=2.;original, it is 2. Tpf index
    
    if not keyword_set(rchisqmax) then rchisqmax = 1000.
    if not keyword_set(spectitle) then spectitle = 'Observed SED of core'
    pwrlaw = keyword_set(rzero)
    if pwrlaw then rz = rzero
    if keyword_set(indsum) then isum = indsum else isum = -1.e35

; Preprocess input data. This includes regridding the images.
    prefit,imagefiles,imagesigma,wavelen,fwhm_beam,fwhm_resolu,back,ncells,racent,deccent, $
	error,SEDFILE=sedfile,XPEAKPIX=xpeakpix,YPEAKPIX=ypeakpix,        $
	ADDNOISE=addnoise,throw_value=throw_value
    if error then return
    infiles = 'regrid_'+imagefiles
    use_sed = keyword_set(sedfile)
    if use_sed then inspec = sedfile
    outlabel = outname

    if keyword_set(addnoise) then begin
    	infiles = 'noisy_'+infiles
    	if use_sed then inspec = 'noisy_'+sedfile
    	label = string(format='("sample",i3.3,"_")',addnoise)
    	outlabel = label+outname
    endif
    
    if keyword_set(makeround) then begin
      make_round,infiles
      infiles = 'round_'+infiles
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
	print,'Writing out '+apmapfile
	writefits,apmapfile,apmap
    endif
    modsed = fltarr(nlambda1,nap)
    obsapflux = fltarr(nbands,nap)
  endif ; if use_sed

; Run a grid of models to get initial parameters for Newton-Raphson step.
; Assume the core is optically thin at the longest wavelength.
    longest = max(lambda1set,nlw)
    ok = where(imageset eq imageset and imageset ne 0, count); original its imageset ne 0, changed by XYC; change back ftom gt 0 by XYC
    if count eq 0 then begin
       message,/CONT,'No valid data'
       return
    endif
    gseed = long(1.e9*(imageset(ok(0))- min(imageset(ok)))/ $
            max(imageset(ok)) - min(imageset(ok)))
    
    if 1 then begin
      ;imageset=imageset>0. ;deleted by XYC, only =0 are ignored
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

while ipass lt maxpasses and notconverged do begin
    count_i=0
    while nlayers/2.^count_i gt 1. do count_i=count_i+1
    if pwrlaw then r0grid = rzero else begin
      r0grid = nlayers/2.^(findgen(count_i+ipass*2)+1.)
    endelse
    rtgrid = nlayers/2.^(findgen(count_i+ipass*2)+1.)
    
    if keyword_set(alpha_preset) then alphagrid=[alpha_preset] else $
      alphagrid = 1.2 + 0.2/precise*findgen(7*precise) + 0.05*ipass*randomn(gseed,7*precise)
    if keyword_set(beta_preset) then beta1grid=[beta_preset] else $
      beta1grid = -2.0 + 0.1/precise*findgen(6*precise) + 0.05*ipass*randomn(gseed,6*precise)
    T0grid = 5. + 5.*findgen(6*precise) + ipass*randomn(gseed,6*precise)
    T1grid = 10. + 5.*findgen(4*precise) + ipass*randomn(gseed,4*precise)
    
    ipass = ipass + 1

    nr = n_elements(r0grid)
    na = n_elements(alphagrid)
    nt = n_elements(T0grid)
    nt1 = n_elements(T1grid)
    nrt = n_elements(rtgrid)
    nb = n_elements(beta1grid)
    
    nstate = 7
    if pwrlaw then nstate = nstate -1		; number of unknowns
    if keyword_set(alpha_preset) then nstate = nstate -1
    if keyword_set(beta_preset) then nstate = nstate -1
    
    if keyword_set(alpha_preset) then alpha_i=1 else alpha_i=0
    if keyword_set(beta_preset) then beta_i=1 else beta_i=0
    
    p = fltarr(nstate)
    phimin = 1.e35
    chimin = 1.e35
    iterations = 0
    finished = 0
    print,'Begin coarse minimization'

    left_counts=nr*nt*nt1*nrt*nb*na

    if (not keyword_set(input_p)) or (ipass gt 1) then begin ;XYC's input p
      for ir=0,nr-1 do for it=0,nt-1 do for it1=0,nt1-1 do for irt = 0,nrt-1 $
    	do for ib=0,nb-1 do begin
          ;if ir eq nr-1 then iamin = na-1 else iamin = 0
          for ia = 0,na-1 do begin;for ia = iamin,na-1 do begin
          	p(1) = beta1grid(ib)
          	p(2-beta_i) = alphagrid(ia)
          	p(3-alpha_i-beta_i) = T0grid(it)
          	p(4-alpha_i-beta_i) = T1grid(it1)
          	p(5-alpha_i-beta_i) = rtgrid(irt)
          	if pwrlaw then r0 = rzero else begin
          	    r0 = r0grid(ir)
          	    p(6-alpha_i-beta_i) = r0
          	endelse
                left_counts=left_counts-1
      ;;if keyword_set(alpha_preset) then p_temp_2=alpha_preset else p_temp_2=alphagrid(ia)
    	if isum gt -1.e30 then gamma1 = isum - alphagrid(ia) else gamma1 = gammainit
    					; index used in temperature profile
      ;nlayers is for the whole sphere, != npixel
    	T = p(4-alpha_i-beta_i) + (p(3-alpha_i-beta_i) - p(4-alpha_i-beta_i))/(1. + (findgen(nlayers)/p(5-alpha_i-beta_i))^gamma1)
    	density = 1./(1. + (findgen(nlayers)/r0)^alphagrid(ia))
    	
    	if max_resolve(nlw) le 0.001 then $
    	  unitmass = coreimage(T,density,1.,rout,a,beta1grid(ib),d, $
                   lambda1set(nlw),nlayers,npixels,npixels) else begin
    	  routset = (findgen(nlayers)+1.)*rout/nlayers
    	  ;rout is in pc = 0.5*npixels*pixel*!dtor*d
    	  density_part=density*(routset le max_resolve(n)*!dtor*d)
    	  T_part=T*(routset le max_resolve(n)*!dtor*d)
    	  unitmass = coreimage(T_part,density_part,1.,rout,a,beta1grid(ib),d, $
                   lambda1set(nlw),nlayers,npixels,npixels)
    	endelse
    	
    	if 0 then p(0) = total(imageset(*,*,nlw))/total(unitmass) else $
    	  p(0) = total(imageset(*,*,nlw))/total(unitmass*(imageset(*,*,nlw) ne 0))
    	chi = phicore(p)
    	if chi lt chimin then begin
    	    pinit = p
    	    print, 'p(now the best) = ', p
            print, 'left counts = ', left_counts
    	    chimin = chi
    	endif
          endfor
        endfor
    
        print,'Initial estimate of core center temperature =',pinit(3-alpha_i-beta_i),' K'
        print,'Initial estimate of core mass        =',pinit(0),' solar masses'
        if keyword_set(beta_preset) then p_temp_1=beta_preset else p_temp_1=pinit(1)
        print,'Initial estimate of beta		=',p_temp_1
    
        p = pinit
    endif ;test p value
    if keyword_set(input_p) and ipass eq 1 then p=input_p ;XYC's input p
    
    phimin = 1.e35
    chi = phicore(p)
    print, 'p(init) = ', p
    printf_p_init=p
    print, 'resuced chi2 of p = ', chi
    print,'Now doing final minimization using Powell procedure'

  printf_p_after_set0=list()
  printf_p_after_set1=list()


  if fit_setting(1) ne 0 then begin;full p
  fit_loop=0
  iterations=0
  while (fit_loop lt fit_setting(1)) or (fit_setting(1) eq -1) do begin
    print,'Now fit full p, loop = ', fit_loop
    iterations = 0
    ftol = 1e-4
    xi = fltarr(nstate,nstate)
    xi(indgen(nstate),indgen(nstate)) = 1.

    if keyword_set(double) then powell,p,xi,ftol,fmin,'phicore',/DOUBLE,ITMAX=1000000 $
        else powell,p,xi,ftol,fmin,'phicore',ITMAX=1000000
    ;print,'xi=',xi
    ;if (fmin le rchisqmax) and (p(4-alpha_i-beta_i) gt 5) then notconverged = 0
    if (fmin le rchisqmax) then notconverged = 0
    print, 'p(final) = ', p
    printf_p_after_set1.add,p
    print,'ipass = ', ipass
    fit_loop=fit_loop+1
    if fit_setting(1) eq -1 then begin
        if iterations le 30 then break
    endif
  endwhile
  endif


  rout0=rout

  if keyword_set(fit_rout) then begin;p_rout
    phicore_quanzhong=1.
    iterations = 0
    ftol = 1e-4
    xi = fltarr(nstate,nstate)
    xi(indgen(nstate),indgen(nstate)) = 1.
    p_no0=fltarr(nstate)
    p_no0(0:nstate-2)=p(1:*)
    p_no0(nstate-1)=rout*2
    rout0=rout
    powell,p_no0,xi,ftol,fmin,'phicore_rout',ITMAX=1000000
    ;print,'xi=',xi
    p(1:*)=p_no0
    ;if (fmin le rchisqmax) and (p(4-alpha_i-beta_i) gt 5) then notconverged = 0
    if (fmin le rchisqmax) then notconverged = 0
    p(0)=p_cloudmass
    rout=p_no0(nstate-1)
    print, 'p(noQZ) = ', p
    print,'ipass = ', ipass
    delvar,phicore_quanzhong
  endif

endwhile

    finished = 1
    chi = phicore(p)
    r = findgen(nlayers)
    
    
    if keyword_set(beta_preset) then p_temp_1=beta_preset else p_temp_1=p(1)
    if keyword_set(beta_preset) then beta_i=1 else beta_i=0
    if keyword_set(alpha_preset) then p_temp_2=alpha_preset else p_temp_2=p(2-beta_i)
    if keyword_set(alpha_preset) then alpha_i=1 else alpha_i=0

    if isum gt -1.e30 then gamma1 = isum - p_temp_2 else gamma1 = gammainit
    ; index used in temperature profile
    Test = p(4-alpha_i-beta_i) + (p(3-alpha_i-beta_i) - p(4-alpha_i-beta_i))/(1. + (r/p(5-alpha_i-beta_i))^gamma1)
    if pwrlaw then r0 = rzero else r0 = p(6-alpha_i-beta_i)
    densest = 1./(1. + (findgen(nlayers)/r0)^p_temp_2)
       
    npeak=coreimage_npeak_idl(Test,densest,p(0),rout,a,p_temp_1,d,lambda1set(nlw),nlayers,npixels,npixels);seehere

;Write out estimated profiles in .csv
    outfile=outlabel + '.csv'
    WRITE_CSV, outfile, [fmin],[p_temp_1],[npeak],[r0*rout/nlayers],[p_temp_2],[p(3-alpha_i-beta_i)],[p(4-alpha_i-beta_i)],[p(5-alpha_i-beta_i)*rout/nlayers], $
      HEADER=['redu_chi_sq','beta','rho0','r0','alpha','T0','T1','rt']

; Write out estimated profiles.
    outfile = outlabel + '.params'
    print,'Writing out '+outfile
    openw,1,outfile
    printf,1,'core mass           =',p(0),' solar masses'
    printf,1,'beta                =',p_temp_1
    printf,1,'alpha               =',p_temp_2
    printf,1,'r0                  =',r0*rout/nlayers,' pc'
    printf,1,'T0                  =',p(3-alpha_i-beta_i),' K'
    printf,1,'T1                  =',p(4-alpha_i-beta_i),' K'
    printf,1,'rt                  =',p(5-alpha_i-beta_i)*rout/nlayers,' pc'
    printf,1,'reduced chi squared =',fmin
    printf,1,'T(rout)             =',p(4-alpha_i-beta_i) + (p(3-alpha_i-beta_i)$
                                  -p(4-alpha_i-beta_i))/(1. + (nlayers/p(5-alpha_i-beta_i))^gamma1)
    printf,1,'Rho(rout)           =',npeak/(1. + (nlayers/r0)^p_temp_2)
    printf,1,' '
    printf,1,'Model density profile (number density of H2):'
    printf,1,' '
    printf,1,'N(r) = N(0)/(1. + (r/r0)^alpha)'
    printf,1,' '
    printf,1,'where N(0) =',npeak,' cm^-3 and r is in units of layers'
    printf,1,'(1 layer =',rout/nlayers,' pc)'
    printf,1,' '
    printf,1,'Model temperature profile:'
    printf,1,' '
    if isum gt -1.e30 then begin
	printf,1,'T(r) = T1 + (T0-T1)/(1. + (r/rt)^gamma)   K'
	printf,1,'    where gamma =',gamma1, $
	';   the sum of density and temperature indices was constrained to be',$
	isum
    endif else printf,1,'T(r) = T1 + (T0-T1)/(1. + (r/rt)^gamma)K, where gamma =',gamma1
    printf,1,' '
    if keyword_set(routset) then printf,1,'input Rout =', routset else $
        printf,1,'automatic Rout =', rout0
    if keyword_set(fit_rout) then printf,1,'calculated Rout =', rout
    printf,1,'nlayers             =',nlayers
    print,rout
    if keyword_set(printf_p_init) then printf,1,'printf_p_init =',printf_p_init
    if keyword_set(printf_p_after_set0) then begin
        for i_p=0,printf_p_after_set0.count()-1 do begin
            printf,1,'P[fit_set(0)] =',printf_p_after_set0[i_p]
        endfor
    endif
    if keyword_set(printf_p_after_set1) then begin
        for i_p=0,printf_p_after_set1.count()-1 do begin
            printf,1,'P[fit_set(1)] =',printf_p_after_set1[i_p]
        endfor
    endif
    printf,1,' '
    printf,1,'Input fit settings:'
    if keyword_set(weight_index) then printf,1,'weight_index  =',weight_index
    if keyword_set(fit_setting) then printf,1,'fit_setting  =',fit_setting
    ;endif else printf,1,'T(r) = T1 + (T0-T1)/(1. + (r/rt)^gamma)   K'
    close,1
    print,'core mass           =',p(0),' solar masses'
    print,'beta1               =',p_temp_1
    print,'alpha               =',p_temp_2
    print,'r0                  =',r0,' layers'
    print,'T0                  =',p(3-alpha_i-beta_i),' K'
    print,'T1                  =',p(4-alpha_i-beta_i),' K'
    print,'rt                  =',p(5-alpha_i-beta_i),' layers'
    print,'peak N(H2)          =',npeak,' cm^-3'
    if isum gt -1.e30 then print,'alpha+gamma         =',isum
    print,'reduced chi squared =',fmin
    print,' '
    ;print,'Estimated radial profiles:'
    ;print,'T:  ',Test
    ;print,'den:',densest
    profiles1 = fltarr(nlayers,2)
    profiles1(*,0) = Test
    profiles1(*,1) = densest
    outfile = outlabel + '_profiles.fits'
    print,'Writing out '+outfile
    mkhdr,hd,profiles1
    sxaddpar,hd,'CDELT2',0.5*npixels*pixel/nlayers
    sxdelpar,hd,'WAVELEN'
    sxaddpar,hd,'COREMASS',p(0)
    sxaddpar,hd,'beta',p_temp_1
    sxaddpar,hd,'RCHISQ',fmin
    writefits,outfile,profiles1,hd

; Write out model images if necessary.
    if keyword_set(modelimages) then outcore,p,outlabel

; Plot SED if necessary.
    if keyword_set(sedplot) then begin
	device = 'ps'
	set_plot, device
	outfile = outlabel + '_sed.ps'
	print,'Writing out file '+outfile
	device, file=outfile
	if not use_sed then nlambda1 = 0
	sed = fltarr(nlambda1+nbands)
	sigsed = fltarr(nlambda1+nbands)
	lambda1 = fltarr(nlambda1+nbands)
	if use_sed then begin
	    maxcover = max(coverage,maxap)
	    for i = 0,nlambda1-1 do begin
		sed(i) = spectra(i,maxap)
		sigsed(i) = sigspectra(i,maxap)
	    endfor
	    lambda1(0:nlambda1-1) = lambda1spec
	    for n = 0,nbands-1 do begin
		sed(nlambda1+n) = obsapflux(n,maxap)
		;print, '################ print sed #################'
		;print, sed
                if n_elements(imagesigma) eq 1 then $
		  sigsed(nlambda1+n) = sed(nlambda1+n)*imagesigma else $
                  sigsed(nlambda1+n) = sed(nlambda1+n)*imagesigma(n)
		lambda1(nlambda1+n) = lambda1set(n)
	    endfor
	endif else begin
	    for n = 0,nbands-1 do begin
		sed(n) = total(imageset(0:nxset(n)-1, 0:nyset(n)-1, n))
		if n_elements(imagesigma) eq 1 then $
                  sigsed(n) = sed(n)*imagesigma else $
                  sigsed(n) = sed(n)*imagesigma(n)
		lambda1(n) = lambda1set(n)
	    endfor
	endelse
	a_1 = findgen(17) * (!pi*2/16.)
	usersym, 0.25*cos(a_1), 0.25*sin(a_1), /fill
	xr = [30., 3000.]
	yr = [min(sed)/3.,max(sed)*5.]
	plot,lambda1,sed,/xlog,/ylog,xrange=xr,/xstyle,yrange=yr, $
	    /ystyle,xtitle='Wavelength [microns]',ytitle='Flux density [Jy]', $
	    title=outlabel,psym=8
	ploterrorbars,lambda1,sed,0*sigsed,sigsed
	npoints = 100
	if use_sed then begin
	    lambda1mod = lambda1(0) + $
	    (lambda1(nlambda1+nbands-1) - lambda1(0))*findgen(npoints)/(npoints-1)
	    modelspectrum,p,lambda1mod,specmod,apspecmod
	    q = where(lambda1mod gt lambda1spec(nlambda1-1), count)
	    lambda1com = fltarr(nlambda1+count)
	    lambda1com(0:nlambda1-1) = lambda1spec
	    lambda1com(nlambda1:nlambda1+count-1) = lambda1mod(q)
	    sedcom = fltarr(nlambda1+count)
	    sedcom(0:nlambda1-1) = modsed(*,maxap)
	    sedcom(nlambda1:nlambda1+count-1) = apspecmod(q,maxap)
	    oplot,lambda1com,sedcom
	endif else begin
	    lambda1mod = 10.^(alog10(xr(0)) + $
		(alog10(xr(1)) - alog10(xr(0)))*findgen(npoints)/(npoints-1))
	    modelspectrum,p,lambda1mod,specmod,apspecmod
	    oplot,lambda1mod,specmod
	endelse
	device,/close

; Output a text file containing model SEDs.
	outfile = outlabel + '_sed.txt'
	print,'Writing out file '+outfile
	openw,20,outfile
	printf,20,'# Spatially integrated SED'
	printf,20,'# ========================'
	printf,20,'#'
	printf,20,'# Square aperture;  width =',npixels*pixel*3600., $
	    ' arcsec  (',2.*rout,' pc)'
	printf,20,'# Center position (RA,Dec)=',racent,deccent,' deg'
	printf,20,'#'
	printf,20,'#        lambda1          flux'
	printf,20,'#       [microns]        [Jy]'
	for n = 0,n_elements(lambda1mod)-1 do printf,20,format='(1x,2f15.4)', $
	    lambda1mod(n),specmod(n)
	printf,20,'#'
	if use_sed then begin
	    printf,20,'#'
	    printf,20,'# SEDs for individual apertures'
	    printf,20,'# ============================='
	    printf,20,'#'
	    printf,20,format= $
		'("# Aperture size parallel to slit      =",f7.2," arcsec")', $
		aperture
	    printf,20,format= $
	    '("# Aperture size perpendicular to slit =",f7.2," arcsec")',slitwid
	    printf,20,'#'
	    printf,20,'# lambda1    aperture fluxes'
	    printf,20,'# [microns]       [Jy]'
	    outform = string(format='("(f8.3,",i2,"f11.5)")',nap)
	    for n = 0,nlambda1-1 do printf,20,format=outform,lambda1spec(n), $
		modsed(n,*)
	    for n = 0,npoints-1 do if lambda1mod(n) gt lambda1spec(nlambda1-1) $
		then printf,20,format=outform,lambda1mod(n),apspecmod(n,*)
	    printf,20,'#'
	    printf,20,'# Observational fluxes within the SED apertures:'
	    for n = 0,nbands-1 do $
		printf,20,format=outform,lambda1set(n),obsapflux(n,*)
	endif
	close,20
    endif

    print,'End COREFIT on '+systime()
    print,' '
    return

    end
