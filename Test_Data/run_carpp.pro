    pro run_carpp

; Run COREFIT on Ori2-2 data and plot the results.

;--------------------------------------------------------------------------
; Set the various required parameter values.
    a = 0.1  ; grain radius [microns]
    d = 414  ; distance [pc]

    ; FITS file names of input images
    imagefiles = ['image-350.fits', 'image-450.fits', 'image-850.fits']

    wavelen = [350., 450., 850.]; observational wavelengths [microns]
    fwhm = [8.5, 8.5, 14.]       ; PSF FWHMs [arcsec]
    imagesigma = 0.05                   ;fractional uncertainty of input images, [0.05, 0.05, 0.05]
    
    back = [0.,0.,0.];back = [0.5, 0.2, 0.05]     ; background levels [Jy/beam]

    xpeakpix_in=[25.0,25.0,25.0] ; the central pixel
    ypeakpix_in=[25.0,25.0,25.0]
    start_pix=1 ; starts from [1,1], as in DS9

    ncells=49                ; the first input image.  The value of NCELLS
                ; should be chosen to include only the core
                ; of interest, i.e. to exclude (as far as
                ; possible) extraneous emission features.
    
    ;routset=0.2 ;pc
    ;sed     = 'input_spectra.dat' ; name of SED file
    outname = 'testdata'            ; root name for the output files
;--------------------------------------------------------------------------

; Run COREFIT.
    carpp, imagefiles,imagesigma,wavelen,fwhm,back,ncells,a,d, $
        outname,sedfile=sed,rchisqmax=rcm,/modelimages,/aperturemap, $
        XPEAKPIX_IN=xpeakpix_in,YPEAKPIX_IN=ypeakpix_in,$
        routset=routset,start_pix=start_pix,input_p=input_p

    end
