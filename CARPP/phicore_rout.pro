    function phicore_rout,p_no0

; Function to be minimized in corefit model.

    common corecom
    common dencom
    common p_for_phicore
    common parameter_range
    
    ;print, 'Print the p values', p_no0
    
    ;cloudmass = p_no0(0)
    r = findgen(nlayers)
    beta1 = p_no0(0)
    
    if keyword_set(alpha_preset) then alpha_i=1 else alpha_i=0
    beta_i=0
    if pwrlaw then r0 = rz else r0 = p_no0(5-alpha_i); > 1.e-10 original

    rout=p_no0(6-alpha_i)
    
    if keyword_set(alpha_preset) then p_temp_2=alpha_preset else p_temp_2=p_no0(1)
    
    if p_temp_2 ne 0. then density = 1./(1. + (r/r0)^p_temp_2) else density = 1.

    if isum gt -1.e30 then gamma1 = isum - p_temp_2 else gamma1 = gammainit
    ;T = p_no0(3-alpha_i) + (p_no0(2-alpha_i)-p_no0(3-alpha_i))/(1. + (r/(p_no0(4-alpha_i)>1.e-10))^gamma1) original
    T = p_no0(3-alpha_i) + (p_no0(2-alpha_i)-p_no0(3-alpha_i))/(1. + (r/(p_no0(4-alpha_i)))^gamma1)
    Tmin = min(T)
    T = T; > 0.1 original
    nbands = n_elements(lambda1set)
    varset1=varset;float(n_elements(varset))/total(varset)*varset
    
        
    if use_sed then begin
	nap = n_elements(sedmask(0,0,*,0))
	nlambda1 = n_elements(sedmask(0,0,0,*))
    endif
    sum = 0.
    wt = 0.
   
    longest = max(lambda1set,nlw)
    unitmass = coreimage(T,density,1.,rout,a,beta1,d, $
    	                     lambda1set(nlw),nlayers,npixels,npixels)
    
    cloudmass = total(imageset(*,*,nlw))/total(unitmass*(imageset(*,*,nlw) ne 0.))
    ;cloudmass = total(imageset(*,*,nlw))/total(unitmass)
    p_cloudmass=cloudmass
    for n = 0,nbands-1 do begin
	model = coreimage(T,density,cloudmass,rout,a,beta1,d,lambda1set(n), $
	    nlayers,npixels,npixels)
	nx = nxset(n)
	ny = nyset(n)
	smod = fltarr(nx,ny)
	smod((nx-1)/2-(npixels-1)/2:(nx-1)/2+(npixels-1)/2,$
	    (ny-1)/2-(npixels-1)/2:(ny-1)/2+(npixels-1)/2)=model ; there was a pixel+1 bug, fixed by XYC
	model = conga(smod,fwhmset(n))
	if use_sed then for iap = 0,nap-1 do obsapflux(n,iap) = $
	    total(imageset(0:nx-1,0:ny-1,n)*sedmask(*,*,iap,nlambda1-1))
	;sum = sum + total((imageset(0:nx-1,0:ny-1,n) - model)^2 / varset1(n)) ;original
	
	;XYC's JiaQuan rchisq
        if not keyword_set(phicore_quanzhong) then qz1=2. else qz1=phicore_quanzhong
	index_array=fltarr(nx,ny)
	cent_pix=(nx-1)/2
	for indx_i = 0, nx-1 do begin
	  for indx_j = 0, ny-1 do begin
	    dist1 = sqrt((indx_i-cent_pix)^2+(indx_j-cent_pix)^2) > 1.
	    if ((dist1 gt cent_pix) or (imageset(indx_i,indx_j,n) eq 0.)) then index_array[indx_i,indx_j]=0 $
	       else index_array[indx_i,indx_j]= 1./dist1^qz1; * (alog10(dist1+0.1)-alog10(dist1))/alog10(1.1);index of the line has no physical meaning
	  endfor
	endfor
	index_array1=float(nx)*float(ny)*index_array/total(index_array); the weird shunxu is to avoid compu error
	sum = sum + total(index_array1 * (imageset(0:nx-1,0:ny-1,n) - model*(index_array ne 0))^2 / varset1(n)) ;XYC's JiaQuan

	wt = wt + total(index_array ne 0);float(nx)*float(ny)
    endfor

  if use_sed then begin ;use_sed begin
    nspec = n_elements(lambda1spec)
    raddist = shift(dist(npixels,npixels),npixels/2,npixels/2)
    out = where(raddist gt npixels/2,nout)
    coverage = fltarr(nap)
    for iap = 0,nap-1 do coverage(iap) = total(sedmask(*,*,iap,0))
    coverage = coverage/max(coverage)

    for n = 0,nspec-1 do begin
	model = coreimage(T,density,cloudmass,rout,a,beta1,d,lambda1spec(n), $
	    nlayers,npixels,npixels)
	for iap = 0,nap-1 do begin
	  if coverage(iap) gt 0 then begin
	    modflux = total(model*sedmask(*,*,iap,n))
	    if coverage(iap) ge 0.5 then begin
		sum = sum + ((modflux - spectra(n,iap))/sigspectra(n,iap))^2
		wt = wt + 1.
	    endif
	    modsed(n,iap) = modflux
	  endif
	endfor
    endfor
  endif;use_sed end

    rchisq = sum/(wt - n_elements(p_no0))
    if finished then return, rchisq

    if keyword_set(printf_p_init) then begin
        if beta1 lt printf_p_init(1)-delta_beta then rchisq = rchisq + ((beta1 - printf_p_init(1))/0.1)^2
        if beta1 gt printf_p_init(1)+delta_beta then rchisq = rchisq + ((beta1 - printf_p_init(1))/0.1)^2
        if p(2-beta_i) lt printf_p_init(2-beta_i)-delta_alpha then rchisq = rchisq + ((p(2-beta_i) - printf_p_init(2-beta_i))/0.1)^2
        if p(2-beta_i) gt printf_p_init(2-beta_i)+delta_alpha then rchisq = rchisq + ((p(2-beta_i) - printf_p_init(2-beta_i))/0.1)^2
	if p(0) lt printf_p_init(0)/delta_cloudmass then rchisq = rchisq + ((p(0) - printf_p_init(0))/0.1)^2
        if p(0) gt printf_p_init(0)*delta_cloudmass then rchisq = rchisq + ((p(0) - printf_p_init(0))/0.1)^2
    endif
    if Tmin lt param_Tmin then rchisq = rchisq + (Tmin/0.01)^2

    if rchisq ne rchisq then rchisq = 1.e10
    if rchisq lt phimin then begin
    	iterations = iterations + 1
    	if iterations mod 10 eq 0 then print,'Iteration',iterations, $
    	    '     reduced chi squared =',rchisq
    	if iterations mod 10 eq 0 then print,'p_no0(now) = ',p_no0
    	phimin = rchisq
    endif

    return, rchisq

    end

