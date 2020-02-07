;+
; Compare two different runs of canoak v3.3 by simple time series
;-
pro main, path=path, year=year, zoom=zoom, dummy=dummy, optimise=optimise, $
          profile=profile, nostop=nostop, wiso=wiso, d13c=d13c, season=season, soil=soil, $
          daily=daily, h2osoil=h2osoil, h2oleaf=h2oleaf, time=time, $
          values=values, limit=limit
  compile_opt idl2
  dl = ([':','\','/'])[(where([ 'macOS', 'windows', 'unix' ] eq strlowcase(!version.os_family)))[0]]
  ;
  if n_elements(path) eq 0 then $
    path = '/Users/matthias/program/canoak/v3.3/model_output/std/' ; Path to the model runs
  path = path + dl ; to be sure
  zoomin = [0,367]              ; zoom into day range, default [0,367]
  if n_elements(zoom) ne 0 then begin
      if n_elements(zoom) eq 2 then zoomin = zoom $
      else print, 'Warning: zoom no 2-element vector: ', zoom
  endif
  timein = [0,24]              ; zoom into time range, default [0,24] h
  if n_elements(time) ne 0 then begin
      if n_elements(time) eq 2 then timein = time $
      else print, 'Warning: time no 2-element vector: ', time
  endif
  ;
  wisoit     = keyword_set(wiso)     ? 1 : 0
  d13cit     = keyword_set(d13c)     ? 1 : 0
  dummyit    = keyword_set(dummy)    ? 1 : 0
  optimiseit = keyword_set(optimise) ? 1 : 0
  profileit  = keyword_set(profile)  ? 1 : 0
  seasonit   = keyword_set(season)   ? 1 : 0
  soilit     = keyword_set(soil)     ? 1 : 0
  dailyit    = keyword_set(daily)    ? 1 : 0
  h2osoilit  = keyword_set(h2osoil)  ? 1 : 0
  h2oleafit  = keyword_set(h2oleaf)  ? 1 : 0
  ;
  valueit    = keyword_set(values)   ? 1 : 0
  limitit    = keyword_set(limit)    ? 1 : 0
  ;
  fixps = 1
  ;
  undef= 9999.
  vsmow = [2005.2e-6, 155.76e-6, 1.]
  o18limit = 1e-9
  o18limit = 2e-9
  ;
  endung = '.csv'
  if n_elements(year) ne 0 then endung = '-' + autostring(fix(year)) + endung
  files1 = ''
  if seasonit then begin
      files1 = [files1, 'season']
      if d13cit then files1 = [files1, 'season_13c']
  endif
  if soilit then files1 = [files1, 'soil']
  if dummyit then files1 = [files1, 'dummy']
  if h2osoilit then files1 = [files1, 'h2osoil']
  if optimiseit then files1 = [files1, 'optimise']
  files1 = files1 + endung
  nfiles1 = n_elements(files1)
  if nfiles1 eq 1 then $
    nfiles1 = 0 $
  else begin
      files1 = files1[1:*]
      nfiles1 = n_elements(files1)
  endelse
  ;
  files2 = ''
  if dailyit then begin
      files2 = [files2, 'daily_ave']
      if d13c then files2 = [files2, 'daily_ave_13c']
  endif
  files2 = files2 + endung
  nfiles2 = n_elements(files2)
  if nfiles2 eq 1 then $
    nfiles2 = 0 $
  else begin
      files2 = files2[1:*]
      nfiles2 = n_elements(files2)
  endelse
  ;
  files3 = ['profile_air', $
            'profile_fluxes']
  if d13cit then files3 = [files3, 'profile_air_13c']
  if wisoit then files3 = [files3, $
                          'profile_air_wiso', $
                          'profile_fluxes_wiso']
  files3 = files3 + endung
  nfiles3 = n_elements(files3)
  ;
  soilisofile = 'h2osoil_wiso'
  soilisofile = soilisofile + endung
  leafisofile = 'h2oleaf_wiso'
  leafisofile = leafisofile + endung
  ;
  path_base = file_basename(path)
  ;
  ; Read data and plot to files
  charsize = 2.0
  xythick = 5
  lthick = 1
  symsize = 0.5
  maxval = abs(undef)-1.
  minval = -1.*abs(undef)+1.
  print, 'Plot run '+path_base
  for i=0, nfiles1-1 do begin
      if not file_test(path+files1[i]) then goto, next_files1
      if file_test(path+files1[i], /zero_length) then goto, next_files1
      print,'  File: ', files1[i]
      ;mc_fread, file=path+files1[i], delimiter=',', skip=1, header=header, data
      mc_fread, file=path+files1[i], skip=1, header=header, data
      ii = where(data eq undef, count) & if count gt 0 then data[ii] = undef
      sdate = long(reform(data[0,*]))
      sday = sdate/10000l
      shr = (sdate-sday*10000l)/100l
      smn = sdate-sday*10000l-shr*100l
      smn = (smn*60l)/100l
      date = double(sday) + double(shr)/24d0 + double(smn)/1440d0
      date = date - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
      date = date + (!dprecision.eps*abs(date) > !dprecision.eps) ; add small amount for re-calc ascii date
      dhr = double(shr) + double(smn-30)/24d0 ; shift half an hour forward because results are mean of hour before
      ii = where(date ge zoomin[0] and date le zoomin[1] and $
                 dhr ge timein[0] and dhr le timein[1]) & date = date[ii] & data = data[*,ii] & sdate = sdate[ii]
      s = size(data)
      if valueit then begin
          print, 'Date'
          print, '   ', autostring(sdate)
      endif
      tmp = strsplit(files1[i],'.',/extract)
      psfile = 'canoak_'+tmp[0]+'.ps'
      mc_opendev, dev='ps', /helvetica, /portrait, file=psfile
      !P.MULTI=[0,1,4]
      mc_symbol, sym=2, fill=0
      for j=1, s[1]-1 do begin
          plot, date, data[j,*], title='!4', xtitle='!4', ytitle='!4'+header[j], $
                charsize=charsize, font=0, min_value=minval, max_value=maxval, $
                xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                psym=8, symsize=symsize, /ynozero
          if valueit then begin
              print, header[j]
              print, '   ', autostring(data[j,*],3)
          endif
      endfor
      mc_closedev
      if fixps then spawn, 'fixps '+psfile, junk, /stderr
      next_files1:
  endfor
  ;
  for i=0, nfiles2-1 do begin
      if not file_test(path+files2[i]) then goto, next_files2
      if file_test(path+files2[i], /zero_length) then goto, next_files2
      print,'  File: ', files2[i]
      ;mc_fread, file=path+files2[i], delimiter=',', skip=1, header=header, data
      mc_fread, file=path+files2[i], skip=1, header=header, data
      ii = where(data eq undef, count) & if count gt 0 then data[ii] = undef
      s = size(data)
      date = long(reform(data[0,*]))
      ii = where(date ge zoomin[0] and date le zoomin[1]) & date = date[ii] & data = data[*,ii]
      tmp = strsplit(files2[i],'.',/extract)
      psfile = 'canoak_'+tmp[0]+'.ps'
      mc_opendev, dev='ps', /helvetica, /portrait, file=psfile
      !P.MULTI=[0,1,4]
      mc_symbol, sym=2, fill=0
      for j=1, s[1]-1 do begin
          plot, date, data[j,*], title='!4', xtitle='!4', ytitle='!4'+header[j], $
                charsize=charsize, font=0, min_value=minval, max_value=maxval, $
                xthick=xythick, ythick=xythick, thick=lthick, line=0,  $
                psym=8, symsize=symsize, /ynozero
          if valueit then begin
              print, header[j]
              print, '   ', autostring(data[j,*],3)
          endif
      endfor
      mc_closedev
      if fixps then spawn, 'fixps '+psfile, junk, /stderr
      next_files2:
  endfor
  ;
  if profileit then begin
      charsize_save = charsize
      xythick_save = xythick
      charsize = 1.
      xythick = 5
      for i=0, nfiles3-1 do begin
          if not file_test(path+files3[i]) then goto, next_files3
          if file_test(path+files3[i], /zero_length) then goto, next_files3
          print,'  File: ', files3[i]
          ;mc_fread, file=path+files3[i], delimiter=',', skip=1, header=header, data
          mc_fread, file=path+files3[i], skip=1, header=header, data
          ii = where(data eq undef, count) & if count gt 0 then data[ii] = undef
          sdate = long(reform(data[0,*]))
          sday = sdate/10000l
          shr = (sdate-sday*10000l)/100l
          smn = sdate-sday*10000l-shr*100l
          smn = (smn*60l)/100l
          date = double(sday) + double(shr)/24d0 + double(smn)/1440d0
          date = date - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
          date = date + (!dprecision.eps*abs(date) > !dprecision.eps) ; add small amount for re-calc ascii date
          dhr = double(shr) + double(smn-30)/24d0 ; shift half an hour forward because results are mean of hour before
          level = fix(reform(data[1,*]))
          ii = where(date ge zoomin[0] and date le zoomin[1] and $
                     dhr ge timein[0] and dhr le timein[1])
          date = date[ii] & data = data[*,ii] & level = level[ii] & sdate = sdate[ii]
          s = size(data)
          if valueit then begin
              print, 'Date'
              print, '   ', autostring(mc_uniq(sdate))
          endif
          alldays = mc_uniq(long(date))
          ndays = n_elements(alldays)
          tmp = strsplit(files3[i],'.',/extract)
          psfile = 'canoak_'+tmp[0]+'.ps'
          mc_opendev, dev='ps', /helvetica, /portrait, file=psfile
          !P.MULTI=[0,2,2]
          ii = where(date eq date[0], nlevel)
          level = indgen(nlevel)+1
          for j=2, s[1]-1 do begin
              data1 = reform(data[j,*])
              for k=0l, ndays-1 do begin                  
                  ii = where(date ge alldays[k] and date le alldays[k]+1)
                  data2 = data1[ii]
                  date2 = date[ii]
                  plot, [0], [0], title='!4'+autostring(fix(alldays[k])), $
                        xtitle='!4'+header[j], ytitle='!4Model level', $
                        charsize=charsize, font=0, yrange=[0, nlevel], /ystyle, $
                        xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                        xrange=[min(data1),max(data1)], /ynozero
                  if valueit then print, '   ', autostring(fix(alldays[k]))
                  hr = mc_uniq(date2)
                  nhr = n_elements(hr)
                  for l=0, nhr-1 do begin
                      jj = where(date2 eq hr[l], count)
                      oplot, data2[jj], level, line=0, thick=0.35*l
                      if valueit then print, '   ', autostring(data2[jj],3)
                  endfor
              endfor
          endfor
          mc_closedev
          if fixps then spawn, 'fixps '+psfile, junk, /stderr
          next_files3:
      endfor
      charsize = charsize_save
      xythick  = xythick_save
  endif
  ;
  if wisoit then begin
      if h2oleafit then begin
          file = leafisofile
          if file_test(path+file) then begin
              if (1-file_test(path+file, /zero_length)) then begin
                  print,'  File: ', file
                  ;mc_fread, file=path+file, delimiter=',', skip=1, header=header, data
                  mc_fread, file=path+file, skip=1, header=header, data
                  header = strtrim(header,2)
                  ii = where(data eq undef, count) & if count gt 0 then data[ii] = undef
                  sdate = long(reform(data[0,*]))
                  sday = sdate/10000l
                  shr = (sdate-sday*10000l)/100l
                  smn = sdate-sday*10000l-shr*100l
                  smn = (smn*60l)/100l
                  date = double(sday) + double(shr)/24d0 + double(smn)/1440d0
                  date = date - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
                  date = date + (!dprecision.eps*abs(date) > !dprecision.eps) ; add small amount for re-calc ascii date
                  dhr = double(shr) + double(smn-30)/24d0 ; shift half an hour forward because results are mean of hour before
                  ii = where(date ge zoomin[0] and date le zoomin[1] and $
                             dhr ge timein[0] and dhr le timein[1]) & date = date[ii] & data = data[*,ii] & sdate = sdate[ii]
                  s = size(data)
                  ; calc delta
                  mc = strmid(header,0,/reverse_offset)
                  mcii = where(mc eq '2', mccount)
                  if mccount eq 0 then message,'Unknown error 1.'
                  data1 = fltarr(s[1]+3*mccount,s[2])
                  data1[0:s[1]-1,0:s[2]-1] = data
                  header1 = strarr(s[1]+3*mccount)
                  header1[0:s[1]-1] = header
                  for i=0, mccount-1 do begin
                      for k=0, 2 do begin
                          header1[s[1]+3*i+k] = 'd'+header[mcii[i]+k]
                          dat = data[mcii[i]+k,*]
                          ddat = (dat/vsmow[k]-1.)*1000.
                          ;if limitit then lim = o18limit*(vsmow[k]/vsmow[0]) else lim = 0.
                          ;lii = where(abs(dat) lt lim, lcount)
                          ;if lcount gt 0 then ddat[lii] = -1000.
                          data1[s[1]+3*i+k,*] = ddat
                      endfor
                  endfor
                  ; plot
                  if valueit then begin
                      print, 'Date'
                      print, '   ', autostring(sdate)
                  endif
                  tmp = strsplit(file,'.',/extract)
                  psfile = 'canoak_'+tmp[0]+'.ps'
                  mc_opendev, dev='ps', /helvetica, /portrait, file=psfile
                  !P.MULTI=[0,1,3]
                  mc_symbol, sym=2, fill=0
                  for j=1, s[1]-1 do begin
                      plot, date, data1[j,*], title='!4', xtitle='!4', ytitle='!4'+header1[j], $
                            charsize=charsize, font=0, min_value=minval, max_value=maxval, $
                            xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                            psym=8, symsize=symsize, /ynozero
                      if valueit then begin
                          print, header1[j]
                          print, '   ', autostring(data1[j,*],3)
                      endif
                  endfor
                  minval1 = max([-1.*abs(undef)+1.,-1000.+1])
                  !P.MULTI=[0,1,3]
                  mc_symbol, sym=2, fill=0
                  for j=s[1], s[1]+3*mccount-1 do begin
                      plot, date, data1[j,*], title='!4', xtitle='!4', ytitle='!4'+header1[j], $
                            charsize=charsize, font=0, min_value=minval1, max_value=maxval, $
                            xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                            psym=8, symsize=symsize, /ynozero
                      if valueit then begin
                          print, header1[j]
                          print, '   ', autostring(data1[j,*],3)
                      endif
                  endfor
                  mc_closedev
                  if fixps then spawn, 'fixps '+psfile, junk, /stderr
              endif
          endif
      endif
      ;
      if h2osoilit then begin
          file = soilisofile
          if file_test(path+file) then begin
              if (1-file_test(path+file, /zero_length)) then begin
                  print,'  File: ', file
                  ;mc_fread, file=path+file, delimiter=',', skip=1, header=header, data
                  mc_fread, file=path+file, skip=1, header=header, data
                  header = strtrim(header,2)
                  ii = where(data eq undef, count) & if count gt 0 then data[ii] = undef
                  sdate = long(reform(data[0,*]))
                  sday = sdate/10000l
                  shr = (sdate-sday*10000l)/100l
                  smn = sdate-sday*10000l-shr*100l
                  smn = (smn*60l)/100l
                  date = double(sday) + double(shr)/24d0 + double(smn)/1440d0
                  date = date - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
                  date = date + (!dprecision.eps*abs(date) > !dprecision.eps) ; add small amount for re-calc ascii date
                  dhr = double(shr) + double(smn-30)/24d0 ; shift half an hour forward because results are mean of hour before
                  ii = where(date ge zoomin[0] and date le zoomin[1] and $
                             dhr ge timein[0] and dhr le timein[1]) & date = date[ii] & data = data[*,ii] & sdate = sdate[ii]
                  s = size(data)
                                ; calc delta
                                ;Time,h2o_isotopes.lost1,h2o_isotopes.lost2,h2o_isotopes.lost3,h2o_isotopes.lost4,flux.transp1,flux.transp2,flux.transp3,flux.transp4,flux.litterevap1,flux.litterevap2,flux.litterevap3,flux.litterevap4,flux.soilevap1,flux.soilevap2,flux.soilevap3,flux.soilevap4,prof.throughfall1,prof.throughfall2,prof.throughfall3,prof.throughfall4,flux.surfrun1,flux.surfrun2,flux.surfrun3,flux.surfrun4,soil.qdrai1,soil.qdrai2,soil.qdrai3,soil.qdrai4,soil.theta01,soil.theta11,soil.theta21,soil.theta31,soil.theta41,soil.theta51,soil.theta61,soil.theta71,soil.theta81,soil.theta91,soil.theta02,soil.theta12,soil.theta22,soil.theta32,soil.theta42,soil.theta52,soil.theta62,soil.theta72,soil.theta82,soil.theta92,soil.theta03,soil.theta13,soil.theta23,soil.theta33,soil.theta43,soil.theta53,soil.theta63,soil.theta73,soil.theta83,soil.theta93,soil.theta04,soil.theta14,soil.theta24,soil.theta34,soil.theta44,soil.theta54,soil.theta64,soil.theta74,soil.theta84,soil.theta94
                  mc = strmid(header,0,/reverse_offset)
                  mcii = where(mc eq '1', mccount)
                  if mccount eq 0 then message,'Unknown error 1.'
                  data1 = fltarr(s[1]+3*mccount,s[2])
                  data1[0:s[1]-1,0:s[2]-1] = data
                  header1 = strarr(s[1]+3*mccount)
                  header1[0:s[1]-1] = header
                  print,'    extend raw data'
                  for i=0, mccount-1 do begin
                      prec = 0.;1*mean(data[mcii[i],*])
                      for k=0, 2 do begin
                          header1[s[1]+3*i+k] = 'd'+header[mcii[i]+k+1]
                          ;if header1[s[1]+3*i+k] eq 'dprof.throughfall3' then stop
                          dat1 = data[mcii[i]+k+1,*]
                          dat2 = data[mcii[i],*]
                          ddat = (mc_division(dat1, dat2, 0., prec=prec)/vsmow[k]-1.)*1000.
                          if limitit then lim = o18limit*(vsmow[k]/vsmow[0]) else lim = 0.
                          lii = where(abs(dat1) lt lim, lcount)
                          if lcount gt 0 then ddat[lii] = -1000.
                          data1[s[1]+3*i+k,*] = ddat
                          ;data1[s[1]+3*i+k,*] = (mc_division(data[mcii[i]+k+1,*], data[mcii[i],*], 0., prec=prec)/vsmow[k]-1.)*1000.
                      endfor
                  endfor
                  ; plot
                  tmp = strsplit(file,'.',/extract)
                  psfile = 'canoak_'+tmp[0]+'.ps'
                  mc_opendev, dev='ps', /helvetica, /portrait, file=psfile
                  print,'    plot raw data'
                  if valueit then begin
                      print, 'Date'
                      print, '   ', autostring(sdate)
                  endif
                  !P.MULTI=[0,1,4]
                  mc_symbol, sym=2, fill=0
                  for j=1, s[1]-1 do begin
                      plot, date, data1[j,*], title='!4', xtitle='!4', ytitle='!4'+header1[j], $
                            charsize=charsize, font=0, min_value=minval, max_value=maxval, $
                            xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                            psym=8, symsize=symsize, /ynozero
                      if valueit then begin
                          print, header1[j]
                          print, '   ', autostring(data1[j,*],3)
                      endif
                  endfor
                  minval1 = max([-1.*abs(undef)+1.,-1000.+1])
                  print,'    plot delta values'
                  !P.MULTI = [0,1,4]
                  mc_symbol, sym=2, fill=0
                  k = 0
                  for j=s[1], s[1]+3*mccount-1 do begin
                      if ((j-s[1]) mod 3) eq 0 then begin
                          k++
                          plot, date, data1[j-s[1]+k,*], title='!4', xtitle='!4', ytitle='!4'+header1[j-s[1]+k], $
                                charsize=charsize, font=0, $
                                min_value=minval, max_value=maxval, $
                                xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                                psym=8, symsize=symsize, /ynozero
                          if valueit then begin
                              print, header1[j-s[1]+k]
                              print, '   ', autostring(data1[j-s[1]+k,*],3)
                          endif
                      endif
                      yy = reform(data1[j,*])
                      iiyy = where(yy ne -1000., count)
                      if count gt 0 then begin
                          plot, date[iiyy], yy[iiyy], title='!4', xtitle='!4', ytitle='!4'+header1[j], $
                                charsize=charsize, font=0, $
                                ;min_value=minval1, max_value=maxval, $
                                xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                                psym=8, symsize=symsize, /ynozero
                      endif else begin
                          plot, [0.], [0.], /nodata, title='!4', xtitle='!4', ytitle='!4'+header1[j], $
                                charsize=charsize, font=0, $
                                ;min_value=minval1, max_value=maxval, $
                                xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                                psym=8, symsize=symsize, /ynozero
                      endelse
                      if valueit then begin
                          print, header1[j]
                          print, '   ', autostring(data1[j,*],3)
                      endif
                  endfor
                  mc_closedev
                  if fixps then spawn, 'fixps '+psfile, junk, /stderr
              endif
          endif
      endif
  endif
  ;
  !P.MULTI=0
  ;
  if not keyword_set(nostop) then stop
  return
end
