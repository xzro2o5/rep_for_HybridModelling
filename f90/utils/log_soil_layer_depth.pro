pro main, skew=skew
  compile_opt idl2
  ;
  nlayers = 10
  depth = 2.0
  skew = 3.
  ;
  a = (findgen(nlayers)+1.)/float(nlayers)
  layers = depth*(exp(skew*a)-1.)/(exp(skew)-1.)
  ;
  plot, a, layers, psym=-1, symsize=1.5
  for i=0, nlayers-1 do begin
      if i le 8 then $
        zahl = '0'+strtrim(string(i+1),2) $
      else $
        zahl = strtrim(string(i+1),2)
      print, 'z_soil'+zahl+' = '+autostring(layers[i],3)+';'
  endfor
  ;
  stop
end
