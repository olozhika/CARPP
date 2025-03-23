    
function plummer_density_gradient, x
  common plummer_cd
  common share
  ;need to define: r_cd,rc,rho_c,alphap
  r_here=sqrt(x^2+r_cd^2)
  plummer_general=rho_c/(1+(r_here/rc)^alphap)
  return, 2*plummer_general
end
