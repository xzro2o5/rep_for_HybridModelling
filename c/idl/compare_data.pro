;+
; Compare runs of canoak v3.3 with data, time series + data vs model
;-
pro main, path=path, year=year, zoom=zoom, $
          timeit=timeit, diffit=diffit, corrit=corrit, nostop=nostop
  compile_opt idl2
  dl = ([':','\','/'])[(where([ 'macOS', 'windows', 'unix' ] eq strlowcase(!version.os_family)))[0]]
  ;
  if n_elements(path) eq 0 then $ ; Path to the model runs
    mpath = '/Users/matthias/program/canoak/v3.3/model_output/std/' $
  else mpath = path
  mpath = mpath + dl ; to be sure
  if n_elements(year) eq 0 then year = 2002
  zoomin = [0,367]              ; zoom into day range, default [0,367]
  if n_elements(zoom) ne 0 then begin
      if n_elements(zoom) eq 2 then zoomin = zoom $
      else print, 'Warning: zoom no 2-element vector: ', zoom
  endif      
  if n_elements(timeit) eq 0 then timeit = 1 ; data & model time series
  if n_elements(diffit) eq 0 then diffit = 0 ; data-model time series
  if n_elements(corrit) eq 0 then corrit = 1 ; data vs model
  ;
  ;
  mundef= 9999.
  dundef= -9999.
  undef = dundef
  fixps = 1
  ;
  ; Pathes to the model runs
  dpath = '/Users/matthias/data/alexander/'
  dfile = 'Hainich_data_'+autostring(fix(year))+'.csv'
  endung = '.13c.18o.csv'
  season_file = 'season' + autostring(fix(year)) + endung
  soil_file = 'soil' + autostring(fix(year)) + endung
  h2osoil_file = 'h2osoil' + autostring(fix(year)) + endung
  daily_file = 'daily_ave' + autostring(fix(year)) + endung
  optimise_file = 'optimise' + autostring(fix(year)) + endung
  ;
  ; Read data and plot to files
  ;   Date/Time,Day,Hour,p..mbar.,T2..degC.,T..degC.,Tpot..K.,Tdew..degC.,rh....,VPmax..mbar.,VPact..mbar.,VPdef..mbar.,sh..g.kg.,H2OC..mmol.mol.,rho..g.m..3.,wv..m.s.,wd..deg.,rain..tot..mm.,Srain..tot..mm.,rain..for..mm.,Srain..for..mm.,Tpyr..degC.,TDR..W.m..2.,TUR..W.m..2.,SWDR..W.m..2.,SWUR..W.m..2.,Albedo....,DDR..W.m..2.,LWDR..W.m..2.,LWUR..W.m..2.,TRAD..degC.,PAR..umol.sm..2.,Rn..W.m..2.,SHF1..W.m..2.,SHF2..W.m..2.,SHF3..W.m..2.,SHF4..W.m..2.,SHF5..W.m..2.,SHFM..W.m..2.,ST002a..degC.,ST005a..degC.,ST015a..degC.,ST030a..degC.,ST050a..degC.,ST002b..degC.,ST005b..degC.,ST015b..degC.,ST030b..degC.,ST050b..degC.,ST002M..degC.,ST005M..degC.,ST015M..degC.,ST030M..degC.,ST050M..degC.,SM01....,SM02....,SM03....,SM04....,SM05....,SM06....,Rn.SHFM,ETP.mm.,ETP.W.m..2.,SM4-6M,soil.mm,NEE,H,LE
  print, 'Run: ', mpath
  print, 'Read data: ', dfile
  mc_fread, file=dpath+dfile, delimiter=',', skip=1, header=header, data
  ii = where(data eq dundef, count)
  if count gt 0 then data[ii] = undef
  ; and take what we want
  ddate = reform(data[0,*])
  dday = reform(data[1,*])
  dhhmm = long(reform(data[2,*]))
  dhr = dhhmm/100
  dmn = dhhmm-dhr*100
  dP = reform(data[3,*])
  ii = where(dP ne undef)
  dP[ii] = dP[ii]*100.           ; Pa
  dt2 = reform(data[4,*])
  dt = reform(data[5,*])
  drh = reform(data[8,*])
  ii = where(drh ne undef)
  drh[ii] = drh[ii]/100.         ; ratio
  dvp = reform(data[10,*])
  dvpd = reform(data[11,*])
  dh2o = reform(data[13,*])
  drain = reform(data[17,*])
  drainsurface = reform(data[19,*])
  dtraddown = reform(data[22,*])
  dtradup = reform(data[23,*])
  dsraddown = reform(data[24,*])
  dsradup = reform(data[25,*])
  ddraddown = reform(data[27,*])
  dlraddown = reform(data[28,*])
  dlradup = reform(data[29,*])
  dalbedo = reform(data[26,*]) ; dsradup/dsraddown
  dpar = reform(data[31,*])
  drnet = reform(data[32,*])
  dgsoil = reform(data[38,*])
  dtsoil2 = reform(data[49,*])
  dtsoil5 = reform(data[50,*])
  dtsoil10 = reform(data[51,*])
  dtsoil30 = reform(data[52,*])
  dtsoil50 = reform(data[53,*])
  dhsoil5 = reform(data[63,*])
  ii = where(dhsoil5 ne undef)
  dhsoil5[ii] = dhsoil5[ii]/100.         ; ratio
  dhsoil15 = reform(data[55,*])
  ii = where(dhsoil15 ne undef)
  dhsoil15[ii] = dhsoil15[ii]/100.         ; ratio
  dhsoil30 = reform(data[56,*])
  ii = where(dhsoil30 ne undef)
  dhsoil30[ii] = dhsoil30[ii]/100.         ; ratio
  dhsoil = reform(data[64,*])
  dnee = reform(data[65,*])
  dh = reform(data[66,*])
  dle = reform(data[67,*])
  ;
  ; Read model and take appropriate
  ;   daytime,netrn,sumrn,sumh,sumle,canps,gpp,transpiration,canresp,soilresp,boleresp,gsoil,ksi,Tleaf,Tlsun,Tlshade,D13,D13_long,CiCa,gs,ave_dair_13C,diff_par,par,ta,u,ustar,rhov,zL,press.Pa,relative.humidity,CO2air,Tair40,CO2air40,Vpd40,Tleaf40,Sun_A40,Sun_gs_mol40,Sun_rbCO240,Sun_cs40,Sun_ci40,Sun_cc40,Sun_wj40,Sun_wc40,Quantum_sun40,Sun_rh40,Sun_vpd40
  print, 'Read model: ', season_file
  mc_fread, file=mpath+season_file, delimiter=',', skip=1, header=header, model
  ii = where(model eq mundef, count)
  if count gt 0 then model[ii] = undef
  ; and take what we want
  mdate = long(reform(model[0,*]))
  mday = mdate/10000l
  mhr = (mdate-mday*10000l)/100l
  mmn = mdate-mday*10000l-mhr*100l
  mdate = double(mday) + double(mhr)/24d0 + double(mmn)/1440d0
  mdate = mdate - 30d0/1440d0 ; shift half an hour forward because results are mean of hour before
  mdate = mdate + (!dprecision.eps*abs(mdate) > !dprecision.eps) ; add small amount for re-calc ascii date
  mrnet = reform(model[1,*])
  mh = reform(model[3,*])
  mle = reform(model[4,*])
  ;mtrans = reform(model[7,*])
  mnee = reform(model[8,*])+reform(model[9,*])+reform(model[10,*])-reform(model[6,*])
  mgsoil = reform(model[11,*])
  mpar = reform(model[22,*])
  mP = reform(model[28,*])
  mrh = reform(model[29,*])
  mvpd = reform(model[33,*])
  ;
  ;   daytime,netrn,soilh,soille,soil.sfc.temp,soil.T.base,soil.T.15cm,flux.transp,flux.evap,prof.throughfall,soil.soil_mm,soil.qinfl,soil.qdrai,soil.gsoil,soil.qtran,soil.surfrun
  print, 'Read model: ', soil_file
  mc_fread, file=mpath+soil_file, delimiter=',', skip=1, header=header, model
  ii = where(model eq mundef, count)
  if count gt 0 then model[ii] = undef
  mhsoil = reform(model[10,*])
  ;   daytime,soil_mm_50
  print, 'Read model: ', optimise_file
  mc_fread, file=mpath+optimise_file, skip=1, header=header, model;, delimiter=','
  ii = where(model eq mundef, count)
  if count gt 0 then model[ii] = undef
  mhsoil = reform(model[1,*])
  ;
  ;   Time,flux.transp1,flux.evap1,prof.throughfall1,flux.surfrun1,soil.qdrai1,soil.theta01,soil.theta11,soil.theta21,soil.theta31,soil.theta41,soil.theta51,soil.theta61,soil.theta71,soil.theta81,soil.theta91,soil.qin01,soil.qin11,soil.qin21,soil.qin31,soil.qin41,soil.qin51,soil.qin61,soil.qin71,soil.qin81,soil.qin91,soil.qout01,soil.qout11,soil.qout21,soil.qout31,soil.qout41,soil.qout51,soil.qout61,soil.qout71,soil.qout81,soil.qout91
  print, 'Read model: ', h2osoil_file
  mc_fread, file=mpath+h2osoil_file, delimiter=',', skip=1, header=header, model
  ii = where(model eq mundef, count)
  if count gt 0 then model[ii] = undef
  mhsoil5 = reform(model[9,*])
  mhsoil15 = reform(model[12,*])
  mhsoil30 = reform(model[14,*])

  comparr = ['rnet','h','le','nee','gsoil','par','P','rh','vpd','hsoil','hsoil5','hsoil15','hsoil30']
  ncomp = n_elements(comparr)

  charsize = 2.0
  xythick = 5
  lthick = 0.25
  symsize = 0.5
  minval = undef+1.

  dir = file_basename(mpath)
  psfile = 'canoak_data_vs_'+dir+'_'+autostring(fix(year))+'.ps'
  print, 'Plot '+psfile
  mc_opendev, dev='ps', /helvetica, /portrait, file=psfile
;  mc_loadct, 34

  ; Time series
  if timeit then begin
      !P.MULTI=[0,1,4]
      mc_symbol, sym=2, fill=0
      for i=0, ncomp-1 do begin
          junk = execute('dat = d'+comparr[i])
          junk = execute('modl = m'+comparr[i])
          ii = where(ddate ge zoomin[0] and ddate le zoomin[1])
          jj = where(mdate ge zoomin[0] and mdate le zoomin[1])
          yrange=[mc_min([dat[ii],modl[jj]], undef=dundef), mc_max([dat[ii],modl[jj]], undef=dundef)]
          ;if (yrange[0] eq dundef) then yrange=[mc_min(modl, undef=undef), mc_max(modl, undef=undef)]
          plot, ddate[ii], dat[ii], title='!4', xtitle='!4', ytitle='!4'+strupcase(comparr[i]), $
                charsize=charsize, font=0, min_value=minval, yrange=yrange, $
                xthick=xythick, ythick=xythick, thick=lthick, line=0, /ynozero
          oplot, mdate, modl, line=1, thick=lthick+0.5, min_value=minval, psym=8, symsize=symsize;, psym=-1
          ;if comparr[i] eq 'le' then $
          ;  oplot, mdate, mtrans, line=2, thick=lthick+0.5, min_value=minval ;, psym=-1
      endfor
  endif
  
  ; Differences
  if diffit then begin
      !P.MULTI=[0,1,4]
      for i=0, ncomp-1 do begin
          junk = execute('dat = d'+comparr[i])
          junk = execute('modl = m'+comparr[i])
          if comparr[i] eq 'le' then ii = where(modl ne undef and modl lt 1e4) else ii = where(modl ne undef)
          jj = where(dat ne undef)
          dat1 = dat[jj]
          date1 = ddate[jj]
          mod1 = interpol(modl[ii], mdate[ii], date1)
          ii = where(date1 ge zoomin[0] and date1 le zoomin[1])
          plot, date1[ii], dat1[ii]-mod1[ii], title='!4', xtitle='!4', ytitle='!4!9D!4'+strupcase(comparr[i]), $
                charsize=charsize, font=0, min_value=minval, $
                xthick=xythick, ythick=xythick, thick=lthick, line=0, /ynozero
      endfor
  endif

  ; model vs data
  if corrit then begin
      !P.MULTI = [0,2,3]
      for i=0, ncomp-1 do begin
          junk = execute('dat = d'+comparr[i])
          junk = execute('modl = m'+comparr[i])
          if comparr[i] eq 'le' then ii = where(modl ne undef and modl lt 1e4) else ii = where(modl ne undef)
          jj = where(dat ne undef)
          dat1 = dat[jj]
          date1 = ddate[jj]
          mod1 = interpol(modl[ii], mdate[ii], date1)
          ii = where(date1 ge zoomin[0] and date1 le zoomin[1])
                                ;res = mc_linfit(dat1[ii], mod1[ii], undef=undef, yfit=yfit)
          res = ladfit(dat1[ii], mod1[ii]) & yfit = res[0] + res[1]*dat1[ii]
          mc_symbol, sym=2, fill=0
          plot, dat1[ii], mod1[ii], title='!4'+strupcase(comparr[i]), xtitle='!4Data', ytitle='!4Model', $
                charsize=charsize, font=0, min_value=minval, $
                xthick=xythick, ythick=xythick, thick=lthick, psym=8, symsize=symsize, $
                xtick_get=xticks, ytick_get=yticks, /ynozero
          oplot, dat1[ii], yfit, line=0, thick=lthick+1., min_value=minval
          xmin = xticks[0] & xmax = xticks[n_elements(xticks)-1]
          ymin = yticks[0] & ymax = yticks[n_elements(yticks)-1]
          xyouts, xmin+0.1*(xmax-xmin), ymax-0.1*(ymax-ymin), /data, $
                  '!4m = '+autostring(res[1],2), charsize=0.6*charsize, font=0
      endfor
  endif
  ;
  mc_closedev
  if fixps then spawn, 'fixps '+psfile, junk, /stderr
  ;
  !P.MULTI=0
  ;
  if not keyword_set(nostop) then stop
  return
end
