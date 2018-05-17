;+
; NAME:
;       kcorrect.pro
; PURPOSE:
;	compute color corrections factors for given array of fluxes,
;       frequencies, uncertainties, & SED function, if 2 components
; CALLING:
;	kcorrect, nu, f, uncf, parfit, perror, modelno=modelno, CCFILES=ccfiles
; INPUTS:
;	nu = array of frequencies in Hz
;	f = array of flux densities in W/m^2/Hz/sr, same length as nu,
;		will be modified in place
;		(yeah, I know, it should be I, not f)
;		(I guess it could be done in other units)
;	uncf = array of flux density uncertainties, same units as f,
;		will also be modified in place
;	parfit = SED fit parameters used to calculate color correction
;	perror = SED fit parameter uncertainties, used to estimate 
;		uncertainty in color correction factor & propagate
;		it to flux density uncertainty
;	modelno =  int, specifies form of MBB used in fitting so
;		additional T components can corrected separately
; KEYWORDS:
;	CCFILES = list of color correction table files; assumes they
;		are also in directory '~/idl/vlibm/filters/' & have
;		the same format as those listed in col 11 of
;		'FilterSpecs.csv' & are sorted by wavelength to match
;		nu (converted to microns) 1-to-1
;               defaults to list in '~/idl/vlibm/filters/FilterSpecs.csv'
; COMMON BLOCKS:
;	none - called from a file that cannot have common blocks
; EXTERNAL CALLS:
;	greybody.pro
;	*Assumes !csi is defined, is in SI units, & contains mass of hydrogen atom
; HISTORY:
;	Written: Rebecca Pitts, 2017.
;---------------------------------------------------

pro kcorrect, nu, f, uncf, parfit, perror, MODELNO=modelno, CCFILES=ccfiles
;; Add parameters and filter spec file names only to FilterSpecs.csv

  wlum = double(1e+6)*!csi.c/nu
  IF n_elements(ccfiles) LE 0 THEN BEGIN
  	t = read_csv('~/idl/vlibm/filters/FilterSpecs.csv',header=colnames,record_start=1)
  	ccfiles = replicate('NA',n_elements(wlum))
	foreach w, wlum, i do ccfiles[i]=t.field11[where(round(t.field04*10)/10.d eq round(w*10)/10.d)]
  ENDIF
  check_order=0
  foreach s, ccfiles, index do check_order+=round(wlum[index]) EQ stregex(s, '[0-9]+', /extract)
  IF total(check_order) EQ 0 THEN BEGIN
     check_order2=0
     foreach s, reverse(ccfiles), index do check_order +=round(wlum[index]) EQ stregex(s, '[0-9]+', /extract)
     IF total(check_order2) EQ 0 THEN MESSAGE, 'ERROR: color correction tables not found or not properly sorted' $
	ELSE ccfiles = reverse(ccfiles)
  ENDIF
  CASE modelno OF
     0:model = 'mbb1opthin'
     1:model = 'mbb1opthick'
     2:model = 'mbb2comp' ;;--> probably not appropriate for Far-IR data
     3:model = 'mbb3comp' ;;--> this one's not working yet, probably insufficient data
     4:model = 'mbb2src' ;;--> 2 components in emission
     ELSE:MESSAGE,"Error: invalid model number" + string(7b),/INFO
  ENDCASE

  IF (modelno GT 1) THEN BEGIN
     IF finite(parfit[6]) THEN BEGIN
	inc = (double(6e12) - double(2e12))/100.d
	nulist = findgen(100.d)*inc + double(2e12)
	fmod = nulist*call_function(model,nulist,parfit)
	tpi = find_local_extrema(nulist,fmod)
	IF tpi EQ !null THEN tp = 50. ELSE $
	IF (n_elements(tpi) GT 1) THEN tp = max(double(1e+6)*double(!csi.c)/nulist[tpi]) ELSE $
	  tp = double(1e+6)*double(!csi.c)/nulist[tpi]
     ENDIF ELSE tp = 0.
  ENDIF ELSE tp = 0.

  FOREACH ccf, ccfiles, index DO BEGIN
     IF (ccf NE 'NA') or stregex(ccf,'[0-9]+',/extract) EQ round(wlum[index]) THEN BEGIN

	print,'reading '+ccf+'...'+strjoin(['~/idl/vlibm/filters',ccf],'/')
	kgrid = read_ascii(strjoin(['~/idl/vlibm/filters',ccf],'/'),comment_symbol='#')
	vec_beta = kgrid.field1[1:-1,0]
	vec_T = kgrid.field1[0,1:-1]
	grid_cc = kgrid.field1[1:-1,1:-1]
        wmin = 2897.77/max(vec_T) ;;min wavelength in um where T & beta can be interpolated 
        IF (modelno GT 2) AND (wlum[index] LT tp) AND (tp LT wmin) THEN BEGIN
           IF (finite(parfit[6]) AND (perror[6] NE 0)) THEN BEGIN 
              Tmp=parfit[6]
              uTmp=perror[6] 
           ENDIF ELSE continue
        ENDIF ELSE BEGIN
           Tmp=parfit[0]
           uTmp=perror[0]
        ENDELSE
        
     IF N_elements(vec_beta) EQ 1 THEN BEGIN
	   spleen = spl_init(reform(vec_T),reform(grid_cc))
	   terpy = [(Tmp-uTmp),Tmp,(Tmp+uTmp)]
	   terp_cc = SPL_INTERP(reform(vec_T), reform(grid_cc), spleen, terpy) ;these are 2nd derivatives
		;; maybe move up to determine which temp comp to use more scientifically
	   IF terp_cc[1] GT 0 THEN f[index] = f[index]/terp_cc[1] $
	     ELSE Message, "ValueError: color correction factor went to 0. Aborting to protect data. Check values of terpy."
	   unc_cc = (terp_cc[1]-terp_cc[0])/2.d ;;liberal error estimate
	   uncf[index] = f[index]*sqrt((uncf[index]/f[index])^2+(unc_cc/terp_cc[1])^2) ;;propagate

	ENDIF ELSE BEGIN
	   gdims = size(grid_cc,/dimensions)
	   grid_beta = cmreplicate(vec_beta,gdims[1])
	   grid_T = reform(transpose(cmreplicate(vec_T,gdims[0])))
	   triangulate, grid_beta, grid_T, trigs, bs
	   terpx = [(parfit[2]-perror[2])-0.1d,(parfit[2]-perror[2]),parfit[2],(parfit[2]+perror[2]),(parfit[2]+perror[2])+0.1d] ;;beta
           terpy = [(Tmp-uTmp)-1d,(Tmp-uTmp),Tmp,(Tmp+uTmp),(Tmp+uTmp)+1d] ;;T

	   ;; trigrid doesn't allow single-point output; all output are in grids of size (len(x) x len(y))
	   ;; len(x) must = len(y)
	   ;; terpx can't all be the same value
	   ;; 1 throwaway x,y pair isn't always enough
	   ;; IDL, why must you hurt me in this way?

	   terp_cc = trigrid(grid_beta,grid_T,grid_cc,trigs,/quintic,xout=terpx,yout=terpy,extrapolate=bs)

           IF terp_cc[2,2] GT 0 THEN f[index] = f[index]/terp_cc[2,2] $
	     ELSE Message, "ValueError: color correction factor went to 0. Aborting to protect data. Check values of terpx, terpy."

           dummy = [terp_cc[1,1],terp_cc[1,3],terp_cc[3,1],terp_cc[3,3]]
	   unc_cc = (max(dummy)-min(dummy))/2.d ;;liberal error estimate
	   uncf[index] = f[index]*sqrt((uncf[index]/f[index])^2+(unc_cc/terp_cc[2,2])^2) ;;propagate
	ENDELSE

     ENDIF ELSE CONTINUE
  ENDFOREACH
END
