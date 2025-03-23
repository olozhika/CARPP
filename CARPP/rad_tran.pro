function rad_tran, nlayer, layers, lambda1
;layers is the Mass,T,density... arrays of each layer, shape[6,nlayer]
;calculate the brightness through a ray with input layers

common share,  a_d, rho_d, pc, q350, h, k, c, msun, jy, pi, mh2, gdr, arcsec, cloud_d
snu_layer = dindgen(nlayer)
for i = 0 , nlayer-1 do begin
    td = layers[0,i]
    beta1 = layers[1, i] ;layers[1, *]=layers[1, anynumber]
    l_layer = layers[2, i]*pc ;layers[2, *]=LOS length
    nh2 = layers[3, i]	;is cal-ed density array (at the middle of the layer)
    q_layer = q_l(q350, lambda1, beta1)
    tau_layer = nh2*mh2/gdr*3.0*q_layer/4.0/a_d/rho_d*l_layer

    if (i eq 0) then snu_layer[i] = bt(lambda1, td)*(1-exp(-1.0*tau_layer)) $
       else snu_layer[i]=snu_layer[i-1]*(exp(-1.0*tau_layer))+bt(lambda1, td)*(1-exp(-1.0*tau_layer))
endfor

;print, 'flux at layer', snu_layer[nlayer-1]
return, snu_layer[nlayer-1];integrate from layer0 to layer nlayer-1
end
