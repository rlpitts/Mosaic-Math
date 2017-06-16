pro pmap_maker, results, HDR=hdr, SAVEPLOT=saveplot, SAVEFITS=savefits, SAVEDIR=savedir, TAIL=tail
    common keys, keylist, par_info
    ;common savecon, savedir, tail

  IF (N_elements(savedir) EQ 0) THEN savedir='~/mosaic/save'
  IF (N_elements(tail) EQ 0) THEN BEGIN
    tail = 'pfix'+strjoin(strcompress(string(transpose(par_info.fixed[*])),/remove_all))
  ENDIF
  IF ( strpos(savedir,'/',/reverse_search) EQ strlen(savedir)-1 ) THEN fname = savedir+tail ELSE $
     fname = savedir+'/'+tail
  imd = size(results.flux,/dimensions)

  IF (N_elements(hdr) GT 0) THEN BEGIN
     cpx = sxpar(hdr,'crpix*')
     cv = sxpar(hdr,'crval*')
     cd11 = sxpar(hdr,'cdelt1') & cd22 = sxpar(hdr,'cdelt2')
     IF (cd11 EQ 0) AND (cd22 EQ 0) THEN BEGIN
        cd11 = sxpar(hdr,'cd1_1') & cd22 = sxpar(hdr,'cd2_2')
     ENDIF
     nax = sxpar(hdr,'naxis*')
     IF total(imd[0:1] EQ nax) EQ 0 THEN BEGIN ;assume symmetrical cropping for now
        cpx[0] = cpx[0] - abs(imd[0]-nax[0])/2.
        cpx[1] = cpx[1] - abs(imd[1]-nax[1])/2.
     ENDIF
     xarr = reverse(((indgen(imd[0]) + 1 - cpx[0]) * abs(cd11) + cv[0] - 360)*(-1))
     ;have to do it in degrees west b/c inverting the x-axis in IDL images is
     ;literally f*cking impossible.
     yarr = (indgen(imd[1]) + 1 - cpx[1]) * abs(cd22) + cv[1]
  ENDIF ELSE BEGIN
     xarr = indgen(imd[0])+1
     yarr = indgen(imd[1])+1
  ENDELSE
  
  labelsetc = transpose([ transpose(keylist),[['Dust Temperature','T','(K)','13'], $
                                              ['Fiducial Frequency','$\nu_0$','(Hz)','13'], $
                                              ['Dust Emissivity Index','$\beta$','','13'], $
                                              ['$H_2$ Column Density','log$N(H_2)$','($m^{-2}$)','17'], $
                                              ['Gas-to-Dust Mass Ratio','$M_{H_2}$/$M_{dust}$','','13'], $
                                              ['Dust Opacity','$\kappa_0$','($m^2$/kg)','13'], $
                                              ['Background Dust Temperature','T','(K)','13'], $
                                              ['Foreground Dust Temperature','T','(K)','13']]])
    ;;^order is [key(save name), plot title, colorbar label, units, RGB_table]

  ;;Do NOT move this inside the loop after it.
  comments = ['Dust Temperature (K)','Fiducial Frequency (Hz)',$
		'Dust Emissivity Index','H_2 Column Density (log Hmol/m^2)',$
		'Gas-to-Dust Mass Ratio','Dust Opacity (m^2/kg)',$
		'2nd Dust Temperature (K)','2nd H_2 Column Density (log Hmol/m^2)']

  parhdr = hdr & punhdr = hdr
  IF keyword_set(savefits) THEN BEGIN ;;start building headers
     IF sxpar(hdr,'naxis3') NE 0 THEN sxdelpar, parhdr, 'naxis3' & sxdelpar, punhdr, 'naxis3'
     FOR ind=0,N_elements(par_info)-1 DO BEGIN
	card = strcompress('PARINIT'+string(ind),/remove_all)
        IF par_info[ind].fixed EQ 0 THEN BEGIN
	   sxaddpar, parhdr, card, 'free param', comments[ind], before='HISTORY'
	   sxaddpar, punhdr, card, 'free param', comments[ind], before='HISTORY'
	ENDIF ELSE BEGIN
	   ;;insert values of fixed parameters
	   sxaddpar, parhdr, card, par_info[ind].value, comments[ind], before='HISTORY'
	   sxaddpar, punhdr, card, par_info[ind].value, comments[ind], before='HISTORY'
	ENDELSE
     ENDFOR
  ENDIF
    ;dof = n_elements(where(par_info[*].fixed EQ 0, /null))
  
  FOR ind=0,N_elements(par_info)-1 DO BEGIN

     IF par_info[ind].fixed EQ 0 THEN BEGIN

        Parr=results.params[*,*,ind]
        errParr=results.perrs[*,*,ind]
        ;;Parr[where3d(results.flux[*,*,-2] LE 0,xind=x,yind=y)]=!values.d_nan
        ;;errParr[where3d(results.flux[*,*,-2] LE 0,xind=x,yind=y)]=!values.d_nan
        rgbtab = long(labelsetc[ind,4])
        
        IF (strpos(savedir,'/',/reverse_search) EQ 0) THEN fpath = savedir+tail+'_'+labelsetc[ind,0]+'.png' $
        ELSE fpath = savedir+'/'+tail+'_'+keylist[ind]+'.png'
        
	maxp = max(Parr,/nan) < par_info[ind].limits[1]
	minp = min(Parr,/nan) > par_info[ind].limits[0]
	Pmap=image(Parr,xarr,yarr,TITLE=labelsetc[ind,1]+' Map', max_value=maxp, min_value=minp, axis_style=2, $
                   layout=[2,1,1], xtitle='GLON West (deg)', ytitle='GLAT (deg)', RGB_TABLE=rgbtab)
        cb = colorbar(TARGET=Pmap,title=labelsetc[ind,2]+' '+labelsetc[ind,3])
	;;emax = max(errParr,/nan) < maxp
        Perrmap=image(errParr,xarr,yarr,TITLE=labelsetc[ind,1]+' Error Map', /current, axis_style=2, $
                      layout=[2,1,2], xtitle='GLON West (deg)', ytitle='GLAT (deg)', RGB_TABLE=rgbtab)
        cb2 = colorbar(TARGET=Perrmap,title='$\sigma$('+labelsetc[ind,2]+') '+labelsetc[ind,3])

	IF keyword_set(saveplot) THEN BEGIN
           Pmap.save,fpath, border=10, resolution=600, /close
           print,"Saved maps as "+fpath
	ENDIF

	IF keyword_set(savefits) THEN BEGIN
	   pfpath = savedir+'/'+tail+'_'+keylist[ind]+'.fits'
	   pfupath = savedir+'/'+tail+'_'+keylist[ind]+'_err.fits'

	   ;; Don't move these QTTY assignments to the hdr building loop.
	   ;; You only want 1 per image
	   sxaddpar, parhdr, 'QTTY', comments[ind], after='naxis2'
	   sxaddpar, punhdr, 'QTTY', 'error in '+comments[ind], after='naxis2'
	   sxaddpar, parhdr, 'ERRMAP', pfupath, after='QTTY'
	   sxaddpar, punhdr, 'PARMAP', pfpath, after='QTTY'

	   writefits, pfpath, Parr, parhdr
	   writefits, pfupath, errParr, punhdr
	   print, 'wrote '+keylist[ind]+' map to '+pfpath
	   print, 'wrote '+keylist[ind]+' error map to '+pfupath
	ENDIF

     ENDIF ELSE IF (ind EQ 1) AND (par_info[1].fixed EQ 1) THEN BEGIN
  	tau0 = (10.d^(results.params[*,*,3])*results.params[*,*,5]/results.params[*,*,4]) * double(2.8 * !const.mH)
	errtau = sqrt( tau0^2 * ( ( alog(10)*results.perrs[*,*,3] )^2 + (results.perrs[*,*,4]/results.params[*,*,4])^2 $
		 + (results.perrs[*,*,5]/results.params[*,*,5])^2 ) )
	fidfreq = strcompress(string(round(par_info[1].value/1e+9)),/remove_all)
	taumap=image(tau0,xarr,yarr,TITLE='Map of Optical Depth at '+fidfreq+' GHz', axis_style=2, $
                   layout=[2,1,1], xtitle='GLON West (deg)', ytitle='GLAT (deg)', RGB_TABLE='13')
	cb = colorbar(TARGET=taumap,title='$\tau_0$')
	tauerrmap=image(errtau,xarr,yarr,TITLE='Error Map of Optical Depth at '+fidfreq+' GHz', axis_style=2, $
			/current, layout=[2,1,2], xtitle='GLON West (deg)', ytitle='GLAT (deg)', RGB_TABLE='13')
	cb = colorbar(TARGET=tauerrmap,title='$\sigma$($\tau_0$)')
     ENDIF
  ENDFOR
END
