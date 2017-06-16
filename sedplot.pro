pro sedplot, results, xpix, ypix, SAVEPLOT=saveplot, TITLE=title, $
	 		XRANGE=xrange, YRANGE=yrange, SPLFIT=splfit

  nu = results.nu
  nufnu = nu*results.flux[xpix,ypix,*]

  IF (tag_exist(results,'errflux') EQ 1) THEN BEGIN
     nuferr = nu*results.errflux[xpix,ypix,*]
     IF nuferr[0] GT nufnu[0] THEN nuferr[0] = nufnu[0]-(1/nu) & yrange = [nufnu[0]/10,max(nufnu)*10]
     p = errorplot(nu,nufnu,nuferr,'Dk3',errorbar_color='r',/sym_filled,/xlog,/ylog)
  ENDIF ELSE p = plot(nu,nufnu,'Dk3',/sym_filled,/xlog,/ylog)

  p.xtitle = '$\nu$ (Hz)'
  p.ytitle = '$\nu$$F_{\nu}$ (W/$m^2$)'
  IF keyword_set(title) THEN p.title = title ELSE $
	p.title = 'SED fit at ('+strcompress(string(xpix))+', '+strcompress(string(ypix))+')'
  IF N_elements(xrange) NE 0 THEN p.xrange = xrange ELSE p.xrange = [1e+11,1e+13]
  IF N_elements(yrange) NE 0 THEN p.yrange = yrange

  nufmod = nu*mbb1opthin(nu,results.params[xpix,ypix,*])

  IF keyword_set(splfit) THEN BEGIN
     inc = ( alog10(max(nu)) - alog10(min(nu)) )/100.
     lognu = findgen(100)*inc + alog10(min(nu))
     logf = alog10(nufmod)
     splf = 10^( spline(alog10(nu),logf,lognu) )
     ps = plot(10^lognu,splf,':2b',/overplot)
  ENDIF

  p1 = plot(nu,nufmod,'X4b',/overplot)

  IF (N_elements(saveplot) GT 0) THEN BEGIN
     IF (isa(saveplot,'string') EQ 1) THEN p.Save, saveplot+'.png', border=10, resolution=500, /CLOSE $
     ELSE p.Save, 'testSED.png', border=10, resolution=500, /CLOSE
     print, 'plot saved to CWD'
  ENDIF
END

