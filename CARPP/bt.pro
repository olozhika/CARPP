function bt, lambda1, td
common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, arcsec, cloud_d

return, 2.0*h*c/lambda1^3.0/( exp(h/k*c/lambda1/td) -1 );ans is in erg/(sr*cm^2)
end
