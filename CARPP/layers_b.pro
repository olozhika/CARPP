pro layers_b, shells, nshell, layers, nlayer, b

;make a layer array for an impact parameter b and known shells
;shells are arranged from inner to outer

if (b gt shells[2, nshell-1]) then begin
	print, 'b is larger than cloud radius!'
	stop
endif


;find the inner most shell that b passes
i_center = 0
while (b gt shells[2,i_center] ) do i_center = i_center+1
nlayer = (nshell-i_center-1)*2+1
layers = dindgen(6, nlayer)
r_center = shells[2, i_center]
;for the back half of the shperes

;print, 'nlayer', nlayer
case nlayer of
1: begin
	layers[*,0] = shells[*,i_center]
	layers[2,0] = 2.0*sqrt(r_center^2-b^2)
	end
else: begin
	;for the back half of the shperes
	for i = 0, (nlayer-1)/2-1 do begin
	i_out = nshell-i-1 ;-1 is bcz it is index i
	i_in = i_out-1
	r_out = shells[2,i_out]
	r_in  = shells[2,i_in]
	layers[*,i] = shells[*,i_out]
	layers[2,i] =sqrt(r_out^2-b^2)-sqrt(r_in^2-b^2) ;Length of sight, which the light passes
	endfor
        ;for the middle shell
	layers[*,(nlayer-1)/2] = shells[*,i_center]
	layers[2,(nlayer-1)/2] = sqrt(r_center^2-b^2)*2.0

	;for the front half of the shperes
	for i=(nlayer-1)/2+1, nlayer-1 do begin
	i_in = i_center+i-((nlayer-1)/2+1)
	i_out = i_in+1
	r_out = shells[2,i_out]
	r_in  = shells[2,i_in]
	layers[*,i] = shells[*,i_out]
	layers[2,i] =sqrt(r_out^2-b^2)-sqrt(r_in^2-b^2)
	endfor	
	end
endcase

end
