    pro modelspectrum,p,lambda1mod,specmod,apspecmod

; Calculate spectrum of corefit model.  Output the spatially integrated 
; spectrum and the spectrum in a set of apertures defined by apmask.

    common corecom

    cloudmass = p(0)
    beta1 = p(1)
    if keyword_set(alpha_preset) then alpha_i=1 else alpha_i=0
    if keyword_set(beta_preset) then beta_i=1 else beta_i=0
    if keyword_set(alpha_preset) then p_temp_2=alpha_preset else p_temp_2=p(2-beta_i)
    if keyword_set(beta_preset) then beta1=beta_preset
    if pwrlaw then r0 = rz else r0 = p(6-alpha_i-beta_i) > 1.e-10
    r = findgen(nlayers)
    if p_temp_2 ne 0. then density = 1./(1. + (r/r0)^p_temp_2) else $
	density = 1.
    if isum gt -1.e30 then gamma1 = isum - p_temp_2 else gamma1 = gammainit
    T = p(4-alpha_i-beta_i) + (p(3-alpha_i-beta_i)-p(4-alpha_i-beta_i))/(1. + (r/p(5-alpha_i-beta_i))^gamma1)
    Tmin = min(T)
    T = T; > 0.1 original
    nbands = n_elements(lambda1mod)
    specmod = fltarr(nbands)

    if use_sed then begin
	nap = n_elements(sedmask(0,0,*,0))
	nlambda1 = n_elements(sedmask(0,0,0,*))
	apspecmod = fltarr(nbands,nap)
    endif

    for n = 0,nbands-1 do begin
      if max_resolve(n) le 0.001 then $
  model = coreimage(T,density,cloudmass,rout,a,beta1,d,lambda1mod(n), $
	    nlayers,npixels,npixels) else begin
	     routset = (findgen(nlayers)+1.)*rout/nlayers
	     ;rout is in pc = 0.5*npixels*pixel*!dtor*d
	     density_part=density*(routset le max_resolve(n)*!dtor*d)
	     T_part=T*(routset le max_resolve(n)*!dtor*d)
	     model = coreimage(T_part,density_part,cloudmass,rout,a,beta1,d,lambda1mod(n), $
	      nlayers,npixels,npixels)
	    endelse
	specmod(n) = total(model)
	if use_sed then for iap = 0,nap-1 do $
		apspecmod(n,iap) = total(model*sedmask(*,*,iap,nlambda1-1))
    endfor

    return

    end
