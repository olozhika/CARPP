
function plummercd, r_cd0,rc0,rho_c0,alphap0,Rout_cm0,bkg0=bkg0 ;all in cm
  common plummer_cd, r_cd, rc, rho_c, alphap, Rout_cm, bkg
  common share

    ; r_cd is a number, all in cm
    r_cd=r_cd0
    rc=rc0
    rho_c=rho_c0
    alphap=alphap0
    Rout_cm=Rout_cm0
    if not keyword_set(bkg0) then bkg=0. else bkg=bkg0
    if r_cd ge Rout_cm then return, bkg $
    else return, QROMO('plummer_density_gradient', 0., sqrt(Rout_cm^2-r_cd^2))+bkg
 end
   