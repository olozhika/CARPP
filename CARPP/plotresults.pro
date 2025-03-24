    pro plotresults

    outname = 'crap'
; Now plot the results.
    coreplot,['Ori1-13_profiles.fits', 'sample001_Ori1-13_profiles.fits',     $
	'sample002_Ori1-13_profiles.fits', 'sample003_Ori1-13_profiles.fits', $
	'sample004_Ori1-13_profiles.fits', 'sample005_Ori1-13_profiles.fits', $
	'sample006_Ori1-13_profiles.fits', 'sample007_Ori1-13_profiles.fits', $
	'sample008_Ori1-13_profiles.fits', 'sample009_Ori1-13_profiles.fits', $
	'sample010_Ori1-13_profiles.fits'], outname+'_results'
    return

    end
