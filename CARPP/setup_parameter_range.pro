pro setup_parameter_range

common parameter_range, delta_beta, delta_alpha, param_Tmin, delta_cloudmass, limit_range, Tindex

Tindex = 2.0 ; T profile index

limit_range = 0 ; set to 1 to limit the parameter ranges

; -----------------------------------------------------------------
; values below are only used when limit_range = 1
; all the delta params can be set to a large number to make it free
delta_beta = 0.1
delta_alpha = 0.2
param_Tmin = 3. ;min(T profile)

delta_cloudmass = 1e9 ; grid/delta_cloudmass < output_cloudmass < delta_cloudmass

end
