pro setup_parameter_range

common parameter_range, delta_beta, delta_alpha, param_Tmin, delta_cloudmass
; all the delta params can be set to a large number to make it free
delta_beta = 0.1
delta_alpha = 0.2
param_Tmin = 3.0

delta_cloudmass = 1e9 ; grid/delta_cloudmass < output_cloudmass < delta_cloudmass

end