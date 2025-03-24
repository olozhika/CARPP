    pro run_carpp

; Run COREFIT on Ori2-2 data and plot the results.

;--------------------------------------------------------------------------
; Set the various required parameter values.
    a = 0.1  ; grain radius [microns]
    d = 414  ; distance [pc]

    ; FITS file names of input images
    imagefiles = ['image-350.fits', 'image-450.fits', 'image-850.fits']

    wavelen = [350., 450., 850.]    ; observational wavelengths [microns]
    fwhm = [8.5, 8.5, 14.]          ; PSF FWHMs [arcsec]
    ;max_resolve_size = [0., 0., 0.005] ; max resolved radius in degree. For interferometric data! 0 means no missing flux in the input image
    imagesigma = [0.05, 0.05, 0.05] ; RMS noise at the core center of input images
    
    back = [0.,0.,0.]               ; background levels [Jy/beam]

    xpeakpix_in=[25.0,25.0,25.0]    ; the central pixel
    ypeakpix_in=[25.0,25.0,25.0]
    start_pix=1                     ; starts from [1,1], as in DS9

    ncells=49    ; used diameter (in pixel) of the first input image
    
    ;routset=0.2 ; full radius of the core [pc]
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
