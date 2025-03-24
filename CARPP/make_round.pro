pro make_round, imagefiles, nsectors=nsectors, grid_pf=grid_pf, percen_list=percen_list
; e.g., make_round, imagefiles, nsectors=24, grid_pf=15, percen_list=[0.8,0.5,0.1,0.05]

;import the fits
;imagefiles = ["regrid_image-450.fits","regrid_image-850.fits"]
;percen_list=[0.7,0.5,0.3,0.1]

ratio_throw=4. ;exceed than throw the direction
grid_percen = 10
grid_use_start=2
grid_use_end=7;grid_percen-1

h = headfits(imagefiles(0))
a = readfits(imagefiles(0),h)
;get fig size & center posi
nx = sxpar(h,'NAXIS1')
ny = sxpar(h,'NAXIS2')

;profile sector & grid settings
if not keyword_set(nsectors) then nsectors = 24
angles = findgen(nsectors) * 360.0 / nsectors

rr_dist_all=findgen(nx,ny,n_elements(imagefiles))
rr_angle_all=findgen(nx,ny,n_elements(imagefiles))
sectored_percen_all = findgen(nsectors,n_elements(imagefiles))
centrl_diffuse=indgen(n_elements(imagefiles))

;main: get ratios
for n=0,n_elements(imagefiles)-1 do begin

  h = headfits(imagefiles(n))
  a = readfits(imagefiles(n),h)
  
  ;get fig size & center posi
  nx = sxpar(h,'NAXIS1')
  ny = sxpar(h,'NAXIS2')
  x_cen = (nx-1)/2 ;starts from 0, not 1, not DS9
  y_cen = (ny-1)/2 ;starts from 0, not 1, not DS9
  
  ;cal dist and angle
  rr_dist=findgen(nx,ny)
  rr_angle=findgen(nx,ny)
  for y_cut=0,ny-1 do begin
    for x_cut=0,nx-1 do begin
      rr_dist[x_cut,y_cut]=sqrt((x_cut-x_cen)^2.0+(y_cut-y_cen)^2.0)
      if rr_dist[x_cut,y_cut] eq 0. then rr_angle[x_cut,y_cut]=114514. else begin
        rr_angle[x_cut,y_cut]=asin((y_cut-y_cen)/rr_dist[x_cut,y_cut])*180.0/!PI
        if x_cut-x_cen lt 0 then rr_angle[x_cut,y_cut]=180.0-rr_angle[x_cut,y_cut]
        rr_angle[x_cut,y_cut]=rr_angle[x_cut,y_cut]+90.0 ; angle starts from 6'o clock
      endelse
    endfor
  endfor
  
  rr_dist_all[*,*,n]=rr_dist
  rr_angle_all[*,*,n]=rr_angle
  
  ;get sectorized pf & percentage
  if not keyword_set(grid_pf) then grid_pf = (nx+1)/2
  xout_temp=findgen(grid_pf+1)/grid_pf
  xout=xout_temp[1:*]*max(rr_dist)/sqrt(2.0)

  if keyword_set(percen_list) then percen=percen_list else begin
    ;increase or decline
    want=where((rr_dist lt max(rr_dist)/1.414) and a ne 0.)
    rr_v=a[want]
    rr_x=rr_dist[want]
    pf=interpol(rr_v,rr_x,xout)

    if pf[round(grid_pf/3)] lt pf[round(grid_pf/2)] then begin ;dark core
      percen_min = max(pf)/a[x_cen,y_cen] ; flux percentage array
      if percen_min le 1. then normal_pf=1
      percen_min=alog10(percen_min)
      percen=(1+(findgen(grid_percen)/(grid_percen-1))*percen_min)[grid_use_start:grid_use_end] ; flux percentage array
      percen=10.^percen
    endif
    
    if (pf[round(grid_pf/3)] ge pf[round(grid_pf/2)]) or keyword_set(normal_pf) then begin ;normal
      percen_min = (pf[-2]/a[x_cen,y_cen])>0.1 ; flux percentage array
      percen_min=alog10(percen_min)
      percen=(0-(findgen(grid_percen)/(grid_percen-1))*(0-percen_min))[grid_use_start:grid_use_end] ; flux percentage array
      percen=10.^percen
      if keyword_set(normal_pf) then delvar,normal_pf
    endif

  endelse
  sectored_percen = findgen(n_elements(percen),nsectors) ; the flux percentage posi
  sectored_pf = findgen(grid_pf,nsectors) ; the radial pf
  for i_cut=0,n_elements(angles)-1 do begin
    angle_want_start=angles[i_cut]
    angle_want_end=angles[i_cut]+angles[1]
    want=where(((rr_angle ge angle_want_start and rr_angle lt angle_want_end) or rr_angle eq 114514.) and a ne 0.)
    rr_v=a[want]
    rr_x=rr_dist[want]
    sectored_percen[*,i_cut]=interpol(rr_x,rr_v,percen*a[x_cen,y_cen]) ; dist_pix
    sectored_pf[*,i_cut]=interpol(rr_v,rr_x,xout) ; calue at the dist
  endfor
  
  ; get upgrid (of dist) ratio
  sectored_percen_mean = findgen(nsectors) ; the flux percentage posi
  for i_cut=0,n_elements(angles)-1 do begin
    sectored_percen_mean[i_cut]=mean(sectored_percen[*,i_cut])
  endfor
  sectored_percen_mean=mean(sectored_percen_mean)/sectored_percen_mean
  sectored_percen_all[*,n]=sectored_percen_mean
  
endfor

; get upgrid (of dist) ratio
use_direction = indgen(nsectors)
sectored_percen_mean = findgen(nsectors)
for n=0,nsectors-1 do begin
  temp_arr=sectored_percen_all[n,*]
  temp_arr=temp_arr[where(temp_arr gt 0.)]
  if max(temp_arr)/min(temp_arr) gt ratio_throw then use_direction[n]=0 $
    else use_direction[n]=1
  sectored_percen_mean[n] = mean(temp_arr)
endfor

;do upgrid
for n=0,n_elements(imagefiles)-1 do begin

  h = headfits(imagefiles(n))
  a = readfits(imagefiles(n),h)

  ;get fig size & center posi
  nx = sxpar(h,'NAXIS1')
  ny = sxpar(h,'NAXIS2')
  x_cen = (nx-1)/2 ;starts from 0, not 1, not DS9
  y_cen = (ny-1)/2 ;starts from 0, not 1, not DS9

  ;cal dist and angle
  rr_dist=rr_dist_all[*,*,n]
  rr_angle=rr_angle_all[*,*,n]

  ;get sectorized pf & percentage
  if not keyword_set(grid_pf) then grid_pf = (nx+1)/2
  xout_temp=findgen(grid_pf+1)/grid_pf
  xout=xout_temp[1:*]*max(rr_dist)/sqrt(2.0)

  if keyword_set(percen_list) then percen=percen_list else begin
    ;increase or decline
    want=where((rr_dist lt max(rr_dist)/1.414) and a ne 0.)
    rr_v=a[want]
    rr_x=rr_dist[want]
    pf=interpol(rr_v,rr_x,xout)

    if pf[round(grid_pf/3)] lt pf[round(grid_pf/2)] then begin ;dark core
      centrl_diffuse[n]=0
      percen_min = max(pf)/a[x_cen,y_cen] ; flux percentage array
      if percen_min le 1. then normal_pf=1
      percen_min=alog10(percen_min)
      percen=(1+(findgen(grid_percen)/(grid_percen-1))*percen_min)[grid_use_start:grid_use_end] ; flux percentage array
      percen=10.^percen
    endif

    if (pf[round(grid_pf/3)] ge pf[round(grid_pf/2)]) or keyword_set(normal_pf) then begin ;normal
      centrl_diffuse[n]=1
      percen_min = (pf[-2]/a[x_cen,y_cen])>0.1 ; flux percentage array
      percen_min=alog10(percen_min)
      percen=(0-(findgen(grid_percen)/(grid_percen-1))*(0-percen_min))[grid_use_start:grid_use_end] ; flux percentage array
      percen=10.^percen
      if keyword_set(normal_pf) then delvar,normal_pf
    endif

  endelse
  sectored_percen = findgen(n_elements(percen),nsectors) ; the flux percentage posi
  sectored_pf = findgen(grid_pf,nsectors) ; the radial pf
  for i_cut=0,n_elements(angles)-1 do begin
    angle_want_start=angles[i_cut]
    angle_want_end=angles[i_cut]+angles[1]
    want=where(((rr_angle ge angle_want_start and rr_angle lt angle_want_end) or rr_angle eq 114514.) and a ne 0.)
    rr_v=a[want]
    rr_x=rr_dist[want]
    sectored_percen[*,i_cut]=interpol(rr_x,rr_v,percen*a[x_cen,y_cen]) ; dist_pix
    sectored_pf[*,i_cut]=interpol(rr_v,rr_x,xout) ; calue at the dist
  endfor

  ; do upgrid (get dist, renew the value to that of the new dist)
  a_new=a
  for y_cut=0,ny-1 do begin
    for x_cut=0,nx-1 do begin
      problem=1
      for i_cut=0,n_elements(angles)-1 do begin
        angle_want_start=angles[i_cut]
        angle_want_end=angles[i_cut]+angles[1]
        if (rr_angle[x_cut,y_cut] ge angle_want_start and rr_angle[x_cut,y_cut] lt angle_want_end) then begin
          problem=0
          break
        endif
      endfor
      if rr_angle[x_cut,y_cut] ne 114514. then begin
        if problem eq 1 then print, 'Why? There might be a problem in ugriding' and return, 0
        if not use_direction[i_cut] then a_new[x_cut,y_cut]=0. else begin
          want_dist=rr_dist[x_cut,y_cut]/sectored_percen_mean[i_cut]
          if want_dist ge max(rr_dist)/1.414 then a_new[x_cut,y_cut]=0. else $
            a_new[x_cut,y_cut]=interpol(sectored_pf[*,i_cut],xout,want_dist)
        endelse
      endif
    endfor
  endfor

  a_new[where(rr_dist ge max(rr_dist)/1.414)]=0.

  ; save a_new to fits file
  ;sxaddpar,h,'UPSCALE',upgrid_ratio ; the full image is upscaled by ...
  outimage='round_'+imagefiles(n)
  writefits,outimage,a_new,h

endfor

; Write out estimated profiles.
outfile = 'make_round.params'
print,'Writing out make_round info'
openw,1,outfile
printf,1,'max permit ratio      =',ratio_throw,' # up/down grid ratio of diff wavelen, directions exceed the value are masked'
printf,1,'is dark core?         =',centrl_diffuse
printf,1,'direction used'
printf,1,use_direction
printf,1,'sectored_percen_mean (larger means scale up)' ;fangda
printf,1,sectored_percen_mean
printf,1,'sectored_percen_all'
printf,1,sectored_percen_all
close,1

end
