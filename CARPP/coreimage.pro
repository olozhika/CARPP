    function coreimage,T,density,cloudmass,rout,a,beta1,d,lambda1,nlayers,nx,ny
;Note 230116: same as 221206 version

; Calculate the spatial intensity distribution of a spherical model cold core.
; Output units: Jy/pixel.
;
; Input parameters:
;	T	=	temperature profile [K], starting from core center
;			(vector array of length NLAYERS)
;	density	=	density profile, relative to core center, i.e.,
;			density(0) = 1.0  (vector array of length NLAYERS)
;	cloudmass=	total mass of cloud [solar masses]
;	rout	=	outer radius of core [pc]
;	a	=	grain radius [microns]
;	beta1	=	index of opacity law (opacity ~ lambda1^beta1)
;	d	=	distance [pc]
;	lambda1	=	wavelength [microns]
;	nlayers	=	number of layers in core model
;	nx,ny	=	dimensions of output image

    common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, $
	arcsec, cloud_d
    common dencom
setup_header,a=a,d=d
    ;a_d = a/1.e4 ;done in setup_header
    ;cloud_d = d*pc ;done in setup_header
    routset = (findgen(nlayers)+1.)*rout/nlayers
    beta1set = fltarr(nlayers) + beta1

; Calculate the number density of gas (neutral hydrogen) molecules as a 
; function of radius.
    ;mH2 = 2.*1.673e-24;already in share			; mass of H_2 molecule [g]
    rcloud = rout*pc	; cloud radius [cm], directly from rout
    rlayer = rcloud/nlayers ;1 layers length in [cm]
    
    if full_sphere then begin ;full sphere
      r = findgen(nlayers)*rlayer/1.e8
      dr = rlayer/1.e8
      NH2 = cloudmass*msun*density/((mh2*1.e24)*(4./3.)*!pi* $
        total(density*((r+dr)^3 - r^3))) ; HN2: /[sum(mass of core)/rho_c]
      npeak = NH2(0)
    endif
    
    ;XYC's integrate column density and rho; modified from XYC's python code "PYandIDL"
    ;the function are in plummer*.pro
    if not full_sphere then begin ; not full sphere
      index_array=fltarr(nx,nx);column density
      cent_pix=(nx-1)/2
      for indx_i = 0, nx-1 do begin
        for indx_j = 0, nx-1 do begin
          dist1 = sqrt((indx_i-cent_pix)^2+(indx_j-cent_pix)^2)
          if (dist1 gt cent_pix) then index_array[indx_i,indx_j]=0. else $
            index_array[indx_i,indx_j]=plummercd_grid(dist1*pc_for_pixel*pc,density,rcloud,rlayer)
        endfor
      endfor
      
      mass_temp=total(index_array)*(pc_for_pixel*pc)^2/Msun*mh2 ;mass in Msun
      
      ;r = findgen(nlayers)*rlayer;in cm
      ;NH2=findgen(nlayers)
      ;for nh2_i=0,nlayers-1 do $
        ;NH2[nh2_i]=plummercd_grid(r[nh2_i],density*cloudmass/mass_temp,rcloud,rlayer)
      NH2=cloudmass*density/mass_temp
      npeak = NH2(0)
    endif
    
; Now calculate core image.
    coregen,T,NH2,cloudmass,routset,beta1set,lambda1,nlayers,nx,cor_image
    return,cor_image
	
    end
