    pro plot_newcorefit

    snrvalues = [100., 50., 20., 10., 5.]
    Nsnr = n_elements(snrvalues)

for i = 0,Nsnr-1 do begin
    snr = snrvalues(i)
    label = string(format='("_",i3.3)',round(snr))
    coreplot,['profiles1','profiles2','profiles3','profiles4','profiles5',  $
	'profiles6']+label+'.fits', '50msun-et-2.0_7wavelen'+label,snr=snr, $
	rchisqmax=1.25

    snr4 = snr*sqrt(7./4.)
    coreplot,['profiles1a','profiles2a','profiles3a','profiles4a',   $
	'profiles5a','profiles6a']+label+'.fits',                    $
	'50msun-et-2.0_4wavelen'+label,snr=snr4,rchisqmax=1.25
endfor

    return

    end
