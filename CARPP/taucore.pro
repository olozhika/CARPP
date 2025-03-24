    function taucore,NH2,a,L

; Calculate normal optical depth of core layer at a reference wavelength
; of 350 microns.
;
; Input parameters:
;	NH2	= molecular hydrogen number density [cm^-3]
;	a	= grain radius [microns]
;	L	= path length [cm]

; Set the values of other parameters.
    mH2 = 2. * 1.673e-24		; mass of H2 molecule [g]
    gtd = 100.				; mass ratio of gas to dust
    rhog = 3.0				; grain density [g cm^-3]
    Q	= 1.36e-4			; absorption efficiency at 350 microns

; Number density of dust:
    Ndust = NH2*mH2/(gtd*(4./3.)*!pi*(1.e-4*a)^3 * rhog)

; Optical depth:
    tau = Q*Ndust*L*!pi*(1.e-4*a)^2

    return,tau

    end
    
