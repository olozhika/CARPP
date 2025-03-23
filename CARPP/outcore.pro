    pro outcore,p,outlabel

; Write out model images for corefit model.

    common corecom

    cloudmass = p(0)
    r = findgen(nlayers)
    beta1 = p(1)
    beta_i=0
    if keyword_set(beta_preset) then begin
      beta_i=1
      beta1=beta_preset
    endif
    alpha_i=0
    p_temp_2=p(2-beta_i)
    if keyword_set(alpha_preset) then begin
      alpha_i=1
      p_temp_2=alpha_preset
    endif
    if pwrlaw then r0 = rz else r0 = p(6-alpha_i-beta_i) > 1.e-10
    if p_temp_2 ne 0. then density = 1./(1. + (r/r0)^p_temp_2) else $
	density = 1.
    if isum gt -1.e30 then gamma1 = isum - p_temp_2 else gamma1 = gammainit
    T = p(4-alpha_i-beta_i) + (p(3-alpha_i-beta_i)-p(4-alpha_i-beta_i))/(1. + (r/(p(5-alpha_i-beta_i)>1.e-10))^gamma1)
    Tmin = min(T)
    T = T; > 0.1
    nbands = n_elements(lambda1set)
    sum = 0.
    wt = 0.

    for n = 0,nbands-1 do begin
	model = coreimage(T,density,cloudmass,rout,a,beta1,d,lambda1set(n), $
	    nlayers,npixels,npixels)
	nx = nxset(n)
	ny = nyset(n)
	smod = fltarr(nx,ny)
	smod((nx-1)/2-(npixels-1)/2:(nx-1)/2+(npixels-1)/2,$
      (ny-1)/2-(npixels-1)/2:(ny-1)/2+(npixels-1)/2)=model ; there was a pixel+1 bug, fixed by XYC
  model = conga(smod,fwhmset(n))
	outfile = outlabel+string(format='("-",i3.3,"_model.fits")',fix(lambda1set(n)))
	print,'Writing out '+outfile
	writefits,outfile,model
    endfor
    return

    end
