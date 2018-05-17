;+
;NAME:
;	sedplot.pro
;PURPOSE:
;	plot SEDs from results of greybody.pro & iter_gb3d.pro
;-
;--------------------------------------

;;helper function:
function wavelen, axis, index, value
   wl = 1e+6*double(!csi.c)/value
   swl= round(wl/floor(alog(wl)))*floor(alog(wl))
   return,string(swl, format='(i8)')
end
;;-------------------------------------
pro sedplot, results, XPIX=xpix, YPIX=ypix, MODELNO=modelno, SAVEPLOT=saveplot, TITLE=title, $
	 		XRANGE=xrange, YRANGE=yrange, LISTP=listp, TEXTPOS=textpos

  nu = double(results.nu)
  dims = size(reform(results.flux),/dimensions)
  IF (N_elements(dims) EQ 1) THEN BEGIN
     xpix=!null & ypix=!null
     nufnu = double(nu*results.flux[*])
     pvec = results.params[*]
     upvec = results.perrs[*]
  ENDIF ELSE BEGIN
     nufnu = double(nu*results.flux[xpix,ypix,*])
     pvec = results.params[xpix,ypix,*]
     upvec = results.perrs[xpix,ypix,*]
  ENDELSE

  CASE modelno OF
    0:model = 'mbb1opthin'
    1:model = 'mbb1opthick'
    2:model = 'mbb2comp' ;;--> prob not appropriate for Far-IR data
    3:model = 'mbb3comp' ;;--> this one's not working yet, probably insufficient data
    4:model = 'mbb2src' ;;--> 2 components in emission
    ELSE:model = 'mbb1opthick' ;;default
  ENDCASE

  IF (modelno GT 2) THEN IF ~finite(pvec[6]) THEN model = 'mbb1opthick'


  IF N_elements(xrange) EQ 0 then xrange = double([3e+11,3e+13])
  IF N_elements(yrange) EQ 0 then yrange = double( [min( nufnu[where(nufnu GT 0)] )/10 , max(nufnu)*10] )
  inc = (xrange[1] - xrange[0] )/100.d
  nulist = findgen(100.d)*inc + xrange[0]
  IF (N_elements(dims) EQ 1) THEN BEGIN
     nufmod = nulist*call_function(model,nulist,results.params[*])
     nufmodp = nu*call_function(model,nu,results.params[*])
  ENDIF ELSE BEGIN
     nufmod = nulist*call_function(model,nulist,results.params[xpix,ypix,*])
     nufmodp = nu*call_function(model,nu,results.params[xpix,ypix,*])
  ENDELSE

  IF (tag_exist(results,'errflux') EQ 1) THEN BEGIN
     IF (N_elements(dims) EQ 1) THEN nuferr = nu*reform(results.errflux) $
	ELSE nuferr = nu*results.errflux[xpix,ypix,*]
     IF nuferr[0] GT nufnu[0] THEN nuferr[0] = nufnu[0]-(1/nu)
     fig = errorplot(nu[where(nufnu GT 0)],nufnu[where(nufnu GT 0)],nuferr[where(nufnu GT 0)], $
		     'Dk4',/sym_filled,errorbar_color='r');,xrange=xrange,yrange=yrange)
  ENDIF ELSE BEGIN
     fig = plot(nu[where(nufnu GT 0)],nufnu[where(nufnu GT 0)],'Dk4',/sym_filled);,xrange = xrange,yrange=yrange)
  ENDELSE

  fig.xtitle = '$\nu$ (Hz)'
  fig.ytitle = '$\nu$$I_{\nu}$ (W m$^{-2}$ sr$^{-1}$)'
  IF keyword_set(title) THEN fig.title = title ELSE BEGIN
     IF (N_elements(xpix) NE 0) AND (N_elements(ypix) NE 0) THEN $
	fig.title = 'SED fit at ('+strcompress(string(xpix))+', '+strcompress(string(ypix))+')' $
     ELSE fig.title = 'Sample SED'
  ENDELSE

  fig.xrange = xrange
  fig.yrange = yrange
  fig.xlog=1
  fig.ylog=1
  fig['axis2'].tickformat='wavelen'
  fig['axis2'].log=1
  fig['axis2'].title='wavelength ($\mu$m)'
  fig['axis2'].showtext=1
  fig['title'].hide=1

  ops = plot(nu,double(nufmodp),'X4b',/overplot)
  ;;print,(nufmodp[where(nufnu GT 0)]-nufnu[where(nufnu GT 0)])/nufnu[where(nufnu GT 0)]
  opm = plot(nulist,nufmod,':2b',/overplot)
  IF ops EQ !null THEN stop

  !except=2
  IF keyword_set(listp) THEN BEGIN
     pdata = string(pvec,format='(D8.1)')
     updata = string(upvec,format='(D8.1)')
     pdata[1] = string(pvec[1]/double(1e+9),format='(D8.1)')
     updata[1] = string(upvec[1]/double(1e+9),format='(D8.1)')
     pdata[3] = string(pvec[3],format='(D8.2)')
     updata[3] = string(upvec[3],format='(D8.2)')
     IF N_elements(pdata) GT 7 THEN BEGIN
	pdata[7] = string(pvec[7],format='(D8.2)')
     	Updata[7] = string(upvec[7],format='(D8.2)')
     ENDIF

     syms =transpose([['T','(K)'], ['$\nu_0$','(GHz)'], ['$\beta$',''], $
			['log N(H$_2$)','(m$^{-2}$)'], ['$M_{H_2}/M_{dust}$',''], $
			['$\kappa_0$','(m$^2$kg$^{-1}$)'],['T$_W$','(K)'],$
			['log N$_w$(H$_2$+H)','(m$^{-2}$)']])

     FOREACH ind, where(upvec GT 0) DO $
	pdata[ind]=strjoin([strtrim(pdata[ind],2),strtrim(updata[ind],2)],'$\pm$',/single)

     symd = size(pdata,/dimensions)
     FOR i=0,symd[0]-1 DO $
	pdata[i]=strtrim(syms[i,0],2)+'='+strtrim(pdata[i],2)+' '+syms[i,1]

     IF N_elements(textpos) EQ 0 THEN textpos=[0.18,0.62]
     ;txt=text(0.18,0.62,strjoin(pdata,'!C'),alignment=0.01,/normal,font_size=11)
     txt=text(textpos[0],textpos[1],strjoin(pdata,'$\n$'),alignment=0.01,/normal,font_size=11)
  ENDIF  

  IF (N_elements(saveplot) GT 0) THEN BEGIN
     IF (isa(saveplot,'string') EQ 1) THEN fig.Save, saveplot+'.png', border=10, resolution=500, /CLOSE $
     ELSE fig.Save, 'testSED.png', border=10, resolution=500, /CLOSE
     print, 'plot saved to CWD'
  ENDIF
END

