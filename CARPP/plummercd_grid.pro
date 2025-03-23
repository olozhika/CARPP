
function plummercd_grid, r_cd,density,rcloud,rlayer
    ;rcloud = cloud radius [cm], directly from rout
    ;rlayer = 1 layers length in [cm], corresponds to nlayer, not pixel
    ;r_cd is a number, all in cm
    if r_cd ge rcloud then return, 0. $
    else begin
      p_nlayers=rcloud/rlayer
      cd_out=0. ; in cm-2
      if r_cd/rlayer-floor(r_cd/rlayer) gt 0.01 then begin
        i=floor(r_cd/rlayer)
        cd_out=cd_out+density[i]*sqrt(((i+1)*rlayer)^2-r_cd^2)
      endif
      for i=ceil(r_cd/rlayer),round(p_nlayers-1) do begin
        cd_out=cd_out+density[i]*(sqrt(((i+1)*rlayer)^2-r_cd^2)-sqrt((i*rlayer)^2-r_cd^2))
      endfor
    endelse
    return, cd_out*2
 end
   
