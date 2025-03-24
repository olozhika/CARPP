	pro coregen, T, density, cloudmass, rout, beta1, lambda1, nshell, dimen, core_img
; Note 230116: same as 221206 version

	;Di Li Apr 30, 2007
	;adopted from model.pro to conform to coreimage input
	;change input layer file to calling parameters
	;Note difference with Ken's code: rout is an array of the outter radii of each layer. Also beta1 is an array
        ;Di Li May 17, 2007: add core_img as one of the parameters and add b_step = rout/(n_b-0.5) in the code

; Calculate the spatial intensity distribution of a spherical model cold core.
; Output units: Jy/pixel.
;
; Input parameters:
;	T	=	temperature profile [K], starting from core center, (vector array of length nshell)
;	density	=	density profile [cc], (vector array of length nshell)
;	cloudmass=	total mass of cloud [solar masses]
;	rout	=	outer radius of core [pc], (vector array of length nshell)

;	beta1	=	index of opacity law (opacity ~ lambda1^beta1), (vector array of length nshell)

;	lambda1	=	wavelength [microns]
;	nshell	=	number of layers in core model
;	dimen	=	dimensions (both x and y assuming nx=ny)of output image, nimage has to be odd, so that the central pixel goes exactly through the spherical center

;Note that the Ken's parameter of a and d are now part of the common_share global variables that should be set once 
	;	d	=	distance [pc] to the cloud, now is cloud_d in cm
	;	a	=	grain radius [microns], now is a_d in cm


common dencom
common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, arcsec, cloud_d

;setup_header

;input the spherical file
;shell structure
;Td   beta1  r nh2 mass tau
;input_layers, infile, nshell, shells, header


;now setup the shell structure from input parameters
shells = dindgen(6, nshell)
shells[0,*] = T
shells[1,*] = beta1
shells[2,*] = rout
shells[3,*] = density ;N(H2)
shells[4,*] = 0
shells[5,*] = 0

rmax = shells[2,nshell-1]; rmax=max(rout)

;set up the b (impact parameter), b is in pc
;only when calculating mass and flux, convert to cm
;bstep is essentially the pixel size. It is set by rout and dimen and
;make the core fit right into the image.

if full_sphere then begin
  n_b = (dimen+1)/2
  bstep = rmax/(n_b-0.5)
  ;note that the last B should be half a stepinside the sphere
  b_array = dindgen(n_b)*bstep ;XYC removed the FOR below, for speed
  ;for i = 0, n_b-1 do begin
  ;b_array[i] = i*bstep
  ;endfor

  ;now calculate brightness at each b
  ;Note the units, b_array in PC, bnu_b in cgs, return fits in Jy/pixel units
  bnu_b = dindgen(n_b)
  lambda1 = lambda1/1.0e4 ; note that in rad_tran, lambda1 is in cm
  for i =0, n_b-1 do begin
    layers_b, shells, nshell, layers, nlayer, b_array[i]
    bnu_b[i] = rad_tran(nlayer, layers, lambda1)
  endfor


  ;generate the image
  center_pix = n_b-1
  core_img = dindgen(dimen, dimen)
  for i = 0, dimen-1 do begin
    for j = 0, dimen-1 do begin
      dist1 = sqrt((i-center_pix)^2+(j-center_pix)^2)*bstep
      if (dist1 gt rmax) then core_img[i,j]=0 else core_img[i,j]= interpol(bnu_b, b_array, dist1)
    endfor
  endfor
  ;convert brightness to flux    Flux = B*omega_pix = B*pixel^2/D^2
  core_img = core_img*bstep*bstep*pc/cloud_d*pc/cloud_d/jy
endif

if not full_sphere then begin
  n_b = (dimen+1)/2; nshell==layers, dimen==npxiels
  bstep = rmax/(n_b-0.5)
  ;note that the last B should be half a stepinside the sphere
  b_array = dindgen(n_b)*bstep

  ;now calculate brightness at each b
  ;Note the units, b_array in PC, bnu_b in cgs, return fits in Jy/pixel units
  bnu_b = dindgen(n_b)
  lambda1 = lambda1/1.0e4 ; note that in rad_tran, lambda1 is in cm
  for i =0, n_b-1 do begin
    layers_b, shells, nshell, layers, nlayer, b_array[i]
    bnu_b[i] = rad_tran(nlayer, layers, lambda1)
  endfor


  ;generate the image
  center_pix = (dimen-1)/2
  core_img = dindgen(dimen, dimen)
  for i = 0, dimen-1 do begin
    for j = 0, dimen-1 do begin
      dist1 = sqrt((i-center_pix)^2+(j-center_pix)^2)*pc_for_pixel
      if (dist1 gt rmax) then core_img[i,j]=0 else core_img[i,j]= interpol(bnu_b, b_array, dist1)
    endfor
  endfor
  ;convert brightness to flux    Flux = B*omega_pix = B*pixel^2/D^2
  core_img = core_img *(pc_for_pixel*pc) /(cloud_d)^2 *(pc_for_pixel*pc) /jy
  ;the weird *(pcxxx)...*(pcxxx) is to avoid ex
endif

end
