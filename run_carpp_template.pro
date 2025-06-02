    pro run_carpp

; Run COREFIT on test data and plot the results.

;--------------------------------------------------------------------------
; Set the various required parameter values.
    a = 0.1  ; grain radius [microns]
    d = 414  ; distance [pc]

    ; FITS file names of input images; the images have to be in [Jy/beam]
    imagefiles = ['image-350.fits', 'image-450.fits', 'image-850.fits']

    wavelen = [350., 450., 850.]    ; observational wavelengths [microns]
    fwhm = [8.5, 8.5, 14.]          ; PSF FWHMs [arcsec]
    ;max_resolve_size = [0., 0., 0.005] ; max resolved radius in degree. Useful for interferometric data! 0 means no missing flux in the input image [NOT TESTED, NOT SURE WHETHER IT WILL WORK]
    imagesigma = [0.05, 0.05, 0.05] ; RMS noise of input images in [noise / max flux of the used data]
    
    back = [0.,0.,0.]               ; background levels [Jy/beam]

    xpeakpix_in=[25.0,25.0,25.0]    ; the central pixel of your core on each imput map
    ypeakpix_in=[25.0,25.0,25.0]
    start_pix=1                     ; set to 1 if your peakpix_in starts from [1,1], as in DS9

    ncells=49    ; used diameter (in pixel) of the first input image
    ;input_p =  [6.17,-1.75,2.00000,10.00000,15.0000,25.0000,6.25000] ;can input your guess to avoid the grid attempts (which will take a long time); but don't use it if your data quality is bad
    ;input_p in  solor_mass, spectral index, density index, Bonnor-Ebert Central flat radius [layers, 50 is the total raidus of the core in most cases], central temperature, environmental temperature, temperture changing radius [layers, 50 is the total raidus of the core in most cases]
    
    ;routset=0.2 ; set full radius of the core [pc], useful when your FOV is not enough
    ;sed     = 'input_spectra.dat'  ; name of your additional SED file
    outname = 'testdata'            ; root name for the output files
;--------------------------------------------------------------------------

; Run COREFIT.
    carpp, imagefiles,imagesigma,wavelen,fwhm,back,ncells,a,d, $
        outname,sedfile=sed,rchisqmax=rcm,/modelimages,/aperturemap, $
        XPEAKPIX_IN=xpeakpix_in,YPEAKPIX_IN=ypeakpix_in,$
        routset=routset,start_pix=start_pix,input_p=input_p, $
        MAX_RESOLVE_SIZE=max_resolve_size

    end
