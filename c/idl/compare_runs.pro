;+
; Compare two different runs of canoak v3.3 by simple time series
;-
pro main, path1=path1, path2=path2, zoom=zoom, year=year, $
          timeit=timeit, diffit=diffit, corrit=corrit, nostop=nostop
  compile_opt idl2
  dl = ([':','\','/'])[(where([ 'macOS', 'windows', 'unix' ] eq strlowcase(!version.os_family)))[0]]
  ;
  if n_elements(path1) eq 0 then $ ; Path to the 1. model run
    path1 = '/Users/matthias/program/canoak/v3.3/model_output/std/'
  path1 = path1 + dl ; to be sure
  if n_elements(path2) eq 0 then $ ; Path to the 2. model run
    path2 = '/Users/matthias/program/canoak/v3.3/model_output/std/'
  path2 = path2 + dl ; to be sure
  if n_elements(year) eq 0 then year = 2002
  zoomin = [0,367]              ; zoom into day range, default [0,367]
  if n_elements(zoom) ne 0 then begin
      if n_elements(zoom) eq 2 then zoomin = zoom $
      else print, 'Warning: zoom no 2-element vector: ', zoom
  endif      
  if n_elements(timeit) eq 0 then timeit = 1 ; data & model time series
  if n_elements(diffit) eq 0 then diffit = 1 ; data-model time series
  if n_elements(corrit) eq 0 then corrit = 1 ; data vs model
  ;
  mundef = 9999.
  undef= -9999.
  fixps = 1
  ;
  endung = '.13c.18o.csv'
  files1 = ['season', $
            'soil', $
            'h2osoil']
  ;files1 = 'h2oiso'
  file2 = 'daily_ave'
  files1 = files1 + autostring(fix(year)) + endung
  file2 = file2 + autostring(fix(year)) + endung
  ;
  nfiles = n_elements(files1)
  path1_base = file_basename(path1)
  path2_base = file_basename(path2)
  ;
  ; Read data and plot to files
  charsize = 2.0
  xythick = 5
  lthick = 1
  symsize = 0.5
  minval = undef+1.
  print, 'Compare '+path1_base+' and '+path2_base
  for i=0, nfiles-1 do begin
      print,'File: ', files1[i]
      mc_fread, file=path1+files1[i], delimiter=',', skip=1, header=header, data1
      ii = where(data1 eq mundef, count) & if count gt 0 then data1[ii] = undef
      s1date = long(reform(data1[0,*]))
      s1day = s1date/10000l
      s1hr = (s1date-s1day*10000l)/100l
      s1mn = s1date-s1day*10000l-s1hr*100l
      date1 = double(s1day) + double(s1hr)/24d0 + double(s1mn)/1440d0
      date1 = date1 - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
      date1 = date1 + (!dprecision.eps*abs(date1) > !dprecision.eps) ; add small amount for re-calc ascii date
      s1 = size(data1)
      ii = where(date1 ge zoomin[0] and date1 le zoomin[1]) & date1 = date1[ii] & data1 = data1[*,ii]
      mc_fread, file=path2+files1[i], delimiter=',', skip=1, header=header, data2
      ii = where(data2 eq mundef, count) & if count gt 0 then data2[ii] = undef
      s2date = long(reform(data1[0,*]))
      s2day = s2date/10000l
      s2hr = (s2date-s2day*10000l)/100l
      s2mn = s2date-s2day*10000l-s2hr*100l
      date2 = double(s2day) + double(s2hr)/24d0 + double(s2mn)/1440d0
      date2 = date2 - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
      date2 = date2 + (!dprecision.eps*abs(date2) > !dprecision.eps) ; add small amount for re-calc ascii date
      ii = where(date2 ge zoomin[0] and date2 le zoomin[1]) & date2 = date2[ii] & data2 = data2[*,ii]
      s2 = size(data2)
      if s1[1] ne s2[1] then message, files1[i]+' have different # of columns in '+path1+' and '+path2
      tmp = strsplit(files1[i],'.',/extract)
      file = 'canoak_'+path1_base+'_vs_'+path2_base+'-'+tmp[0]+'.ps'
      mc_opendev, dev='ps', /helvetica, /portrait, file=file
      if timeit then begin
          print,'  time'
          !P.MULTI=[0,1,4]
          mc_symbol, sym=2, fill=0
          for j=1, s1[1]-1 do begin
              plot, date1, data1[j,*], title='!4', xtitle='!4', ytitle='!4'+header[j], $
                    charsize=charsize, font=0, min_value=minval, $
                    xthick=xythick, ythick=xythick, thick=lthick, line=0, /ynozero
              oplot, date2, data2[j,*], line=1, psym=8, symsize=symsize, thick=lthick+1, min_value=minval
          endfor
      endif
      if diffit then begin
          print,'  difference'
          !P.MULTI=[0,1,4]
          mc_symbol, sym=2, fill=0
          for j=1, s1[1]-1 do begin
              dat1 = reform(data1[j,*])
              dat2 = reform(data2[j,*])
              ii = where(dat1 ne undef and finite(dat1) eq 1, icount)
              jj = where(dat2 ne undef and finite(dat2) eq 1, jcount)
              if icount gt 1 and jcount gt 1 then begin
                  dat = interpol(dat2[jj], date2[jj], date1[ii])
                  plot, date1[ii], dat-dat1[ii], title='!4', xtitle='!4', ytitle='!4!9D!4'+header[j], $
                        charsize=charsize, font=0, min_value=minval, $
                        xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                        psym=8, symsize=symsize, /ynozero
              endif
          endfor
      endif
      if corrit then begin
          print,'  correlation'
          !P.MULTI=[0,2,3]
          mc_symbol, sym=2, fill=0
          for j=1, s1[1]-1 do begin
              dat1 = reform(data1[j,*])
              dat2 = reform(data2[j,*])
              ii = where(dat1 ne undef and finite(dat1) eq 1, icount)
              jj = where(dat2 ne undef and finite(dat2) eq 1, jcount)
              if icount gt 1 and jcount gt 1 then begin
                  dat = interpol(dat2[jj], date2[jj], date1[ii])
                  ; res = mc_linfit(dat1[ii], dat, undef=undef, yfit=yfit)
                  res = ladfit(dat1[ii], dat) & yfit = res[0] + res[1]*dat1[ii]
                  plot, dat1[ii], dat, title='!4'+header[j], xtitle='!4'+path1_base, ytitle='!4'+path2_base, $
                        charsize=charsize, font=0, min_value=minval, $
                        xthick=xythick, ythick=xythick, thick=lthick, psym=8, symsize=symsize, $
                        xtick_get=xticks, ytick_get=yticks, /ynozero
                  oplot, dat1[ii], yfit, line=0, thick=lthick, min_value=minval
                  xmin = xticks[0] & xmax = xticks[n_elements(xticks)-1]
                  ymin = yticks[0] & ymax = yticks[n_elements(yticks)-1]
                  xyouts, xmin+0.1*(xmax-xmin), ymax-0.1*(ymax-ymin), /data, $
                          '!4m = '+autostring(res[1],2), charsize=0.6*charsize, font=0
              endif
          endfor
      endif
      mc_closedev
      if fixps then spawn, 'fixps '+file, junk, /stderr
  endfor
  ;
  print,'File: ', file2
  mc_fread, file=path1+file2, delimiter=',', skip=1, header=header, data1
  ii = where(data1 eq mundef, count) & if count gt 0 then data1[ii] = undef
  s1 = size(data1)
  date1 = long(reform(data1[0,*]))
  mc_fread, file=path2+file2, delimiter=',', skip=1, header=header, data2
  ii = where(data2 eq mundef, count) & if count gt 0 then data2[ii] = undef
  s2 = size(data2)
  date2 = long(reform(data2[0,*]))
  if s1[1] ne s2[1] then message, files1[i]+' have different # of columns in '+path1+' and '+path2
  tmp = strsplit(file2,'.',/extract)
  file = 'canoak_'+path1_base+'_vs_'+path2_base+'-'+tmp[0]+'.ps'
  mc_opendev, dev='ps', /helvetica, /portrait, file=file
  if timeit then begin
      print,'  time'
      !P.MULTI=[0,1,4]
      mc_symbol, sym=2, fill=0
      for j=1, s1[1]-1 do begin
          plot, date1, data1[j,*], title='!4', xtitle='!4', ytitle='!4'+header[j], $
                charsize=charsize, font=0, min_value=minval, $
                xthick=xythick, ythick=xythick, thick=lthick, line=0, /ynozero
          oplot, date2, data2[j,*], line=1, psym=8, symsize=symsize, thick=lthick+3, min_value=minval
      endfor
  endif
  if diffit then begin
      print,'  difference'
      !P.MULTI=[0,1,4]
      mc_symbol, sym=2, fill=0
      for j=1, s1[1]-1 do begin
          dat1 = reform(data1[j,*])
          dat2 = reform(data2[j,*])
          ii = where(dat1 ne undef and finite(dat1) eq 1, icount)
          jj = where(dat2 ne undef and finite(dat2) eq 1, jcount)
          if icount gt 1 and jcount gt 1 then begin
              dat = interpol(dat2[jj], date2[jj], date1[ii])
              plot, date1[ii], dat-dat1[ii], title='!4', xtitle='!4', ytitle='!4!9D!4'+header[j], $
                    charsize=charsize, font=0, min_value=minval, $
                    xthick=xythick, ythick=xythick, thick=lthick, line=0, $
                    psym=8, symsize=symsize, /ynozero
          endif
      endfor
  endif
  if corrit then begin
      print,'  correlation'
      !P.MULTI=[0,2,3]
      mc_symbol, sym=2, fill=0
      for j=1, s1[1]-1 do begin
          dat1 = reform(data1[j,*])
          dat2 = reform(data2[j,*])
          ii = where(dat1 ne undef and finite(dat1) eq 1, icount)
          jj = where(dat2 ne undef and finite(dat2) eq 1, jcount)
          if icount gt 1 and jcount gt 1 then begin
              dat = interpol(dat2[jj], date2[jj], date1[ii])
              ; res = mc_linfit(dat1[ii], dat, undef=undef, yfit=yfit)
              res = ladfit(dat1[ii], dat) & yfit = res[0] + res[1]*dat1[ii]
              plot, dat1[ii], dat, title='!4'+header[j], xtitle='!4'+path1_base, ytitle='!4'+path2_base, $
                    charsize=charsize, font=0, min_value=minval, $
                    xthick=xythick, ythick=xythick, thick=lthick, psym=8, symsize=symsize, $
                    xtick_get=xticks, ytick_get=yticks, /ynozero
              oplot, dat1[ii], yfit, line=0, thick=lthick, min_value=minval
              xmin = xticks[0] & xmax = xticks[n_elements(xticks)-1]
              ymin = yticks[0] & ymax = yticks[n_elements(yticks)-1]
              xyouts, xmin+0.1*(xmax-xmin), ymax-0.1*(ymax-ymin), /data, $
                      '!4m = '+autostring(res[1],2), charsize=0.6*charsize, font=0
          endif
      endfor
  endif
  mc_closedev
  if fixps then spawn, 'fixps '+file, junk, /stderr
  ;
  !P.MULTI=0
  ;
  if not keyword_set(nostop) then stop
  return
end
