pro setup_header,a=a,d=d

common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, arcsec, cloud_d

pc = 3.08568e18 ;cm
a_d = a*1e-4 ;cm
rho_d = 3.0
q350 = 1.36e-4
cloud_d = double(d*pc)
gdr = 100.

h = 6.626e-27
k = 1.38e-16
c = 3e10
msun = 1.989e33 ;g
jy = 1.0e-23
pi = 3.14159265
mh2 = 1.673e-24*2.0
arcsec = cloud_d/206265.; arcsec corresponding to cm at the cloud distance

end
