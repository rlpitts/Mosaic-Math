;+
; NAME:
;	greybody.pro
; PURPOSE:
;	Compute and fit greybody emission spectrum to given flux data
;	at for an array of wavelengths.
;	Note: lots of helper functions ahead of the main one!
;	Uses the following relationships:
;		kappa = kappa_0*(wl_0/wl)^beta
;		kappa_0*N_H = tau_0, where N = col. density & tau << 1
;	If tau = absorption cross-section & N_H = # of H atoms, then
;		kappa*mu*m_H/G2D = tau/N_H = (tau_0/N_H)*(wl_0/wl)^beta
;	     -->tau = N_H*kappa_0*(mu*m_H/G2D)*(nu/nu_0)^beta
;	but that may make N_H unmanageably big, & conversion to those units hard
; CALLING:
;	gb = greybody( x, fdat, ... see kwarg list)
;	***Some kwargs can be a float, or 2-, or 3-tuple:
;	arg = init_value -OR-
;	arg = [init_value, bool_fixed] -OR-
;	arg = [init_value, min, max]
;	      (use -1 for min/max to switch to default min/max respectively)
; INPUTS:
;	x = array of wavelengths or frequencies, in any units allowed for XUNIT.
;		defaults to wavelength in microns
;	fluxarr = 1D array of fluxes, or stack of images, at each wavelength
;		in either Jy/arcsec^2, MJy/sr, Jy/beam, or W/m^2/hz
; KEYWORDS:
;	ERRARR = float array of errors in flux, same units as fluxarr;
;		overrides weight if both are present
;	WEIGHT = float array of weights for each flux, nominally in 
;		units of 1/flux; backup for errarr; defaults to 
;		flux- & wavelength-weighted scheme I came up with
;	XUNIT = string units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = alphanumeric string units of fdat
;		accepts 'MJysr', 'Jyasec2', or 'Wm2hz' 
;		default is 'Jyasec2'
;	THRESHOLD = lower limit of valid flux values
;		default is 1e-05, assuming FUNIT='Jyarcsec2'
;		will incorporate upper limits later
;	MODELNO = int, used to select form of MBB to fit
;	PLOT = invoke sedplot.pro to plot data & resulting fit
;	DISTANCE = float, distance in pc to objects of interest (would like to be able to
;		supply an array for multiple distances but I'm not that advanced yet...)
;	AUTOCALL = Boolean switch to distinguish solo calls from calls by greybody3D.pro
;	FILTERS = list of filters for which to do color correction (not available with autocall);
;		see 
;  ***See gbpstruct.pro for parameter descriptions
; OUTPUTS:
;	RESULT = either an array of best-fit parameters for SED fit, -or-
;		struct of arrays: nu in Hz, F(nu)_data in W/m^2/Hz, best-fit
;		SED parameters & their errors, & F(nu)_model in W/m^2/Hz 
; COMMON BLOCKS:
;	None b/c I'm a Python native & these things bug the scheisse out of me
;	(I kid - it's a neat trick but I'm writing most of these
;	modules to work alone OR together)
; EXTERNAL CALLS:
;	gbpstruct.pro
;	greybody3D.pro
;	*Assumes !csi is defined, is in SI units, & contains mass of hydrogen atom
;---------------------------------------------------

function planck, nu, T
    	hc2 =  2 * double(!csi.h) / double(!csi.c)^2
    	hkt = double(!csi.h) / ( double(!csi.k) * T )
    	bb = hc2 * double(nu)^3 / (exp( hkt * nu ) - 1)
	return, bb
  end
;---------------------------------------------------

;;modelno=0
function mbb1opthin, X, P
	;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, P[5] = kappa_0
	;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
	;;from Arthur Cox 2000 reproduction of "Allen's Astrophysical Quantities" (1955) - everyone uses this
	;;alternative cited by Wikipedia doesn't seem to be reliable
	;; assume mu = 2.8 & omega accounted for by xunit coversion

	Inu = planck(X, P[0]) * ( (X / P[1])^P[2] ) * ( 10.D^P[3] * 2.8 * double(!csi.mH) / P[4] ) * P[5]
	return, Inu
  end
;---------------------------------------------------

;;modelno=1
;; ASSUMPTIONS: negligible background - all emission local to source & self-attenuating
function mbb1opthick, X, P
	;; P[0] = Tc, P[1] = nu_0, P[2] = beta, P[3] = log(N_col_c), P[4] = G2D, P[5] = kappa_0
	;; assume mu = 2.8 & no need for omega
	;; PARAMETERS CANNOT BE REARRANGED (yes, it sucks)

	bb = planck(X, P[0])
	tau = ( (X / P[1])^P[2] ) * ( 10.D^P[3] * 2.8 * double(!csi.mH) / P[4] ) * P[5]
    	Inu = bb * ( 1 - exp(-tau) )
	return, Inu
  end
;---------------------------------------------------


;;modelno=2
;; probably shouldn't use separately; should be equivalent to mbb1opthick(T)-mbb1opthick(T_bg)
function mbb2comp, X, P
        ;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, 
        ;; P[5] = kappa_0, P[6] = T_bg
        ;; (see e.g. Magnum & Shirley 2015 review)
        ;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
        ;; assume mu = 2.8, UNIFORM DUST COMPOSITION* & omega accounted for by xunit coversion
        ;;*This is a whopper if there ever was one

	bb = planck(X, P[0])
	tau = ((X / P[1])^P[2] ) * ( 10.D^P[3] * 2.8 * double(!csi.mH) / P[4] ) * P[5]
	bb_bg = planck(X, P[6])
	Inu = ( bb * (1 - exp(-tau)) ) + ( bb_bg * exp(-tau) )
        return, Inu
  end
;---------------------------------------------------

;;modelno=3
function mbb3comp, X, P ;deriv, dP 
;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, 
;; P[5] = kappa_0, P[6] = T_bg, P[7] = T_fg

;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
;; assume mu = 2.8, UNIFORM DUST COMPOSITION* & omega accounted for by xunit coversion
;;*again, majorly suspect assumption

	bb = planck(X, P[0])
	tau = ((X / P[1])^P[2] ) * ( 10.D^P[3] * 2.8 * double(!csi.mH) / P[4] ) * P[5]
	bb_bg = planck(X, P[6])
	bb_fg = planck(X, P[7])
	Inu = ( bb * (1 - exp(-tau)) ) + ( bb_bg * exp(-tau) ) + bb_fg
	return, Inu
  end
;---------------------------------------------------


;;modelno=4
;; probably shouldn't use separately; should be equivalent to mbb1opthick(T)-mbb1opthick(T_bg)
function mbb2src, X, P
        ;; P[0] = T_c, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, 
        ;; P[5] = kappa_0, P[6] = T_h, P[7] = log(N_col)_h
        ;; (see e.g. Magnum & Shirley 2015 review)
        ;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
        ;; assume mu = 2.8, UNIFORM DUST COMPOSITION* & omega accounted for by xunit coversion
        ;;*This is a whopper if there ever was one

	bbc = planck(X, P[0])
	tauc = ((X / P[1])^P[2] ) * ( 10.D^P[3] * 2.8 * double(!csi.mH) / P[4] ) * P[5]
	bbw = planck(X, P[6])
	tauw = ((X / P[1])^P[2] ) * ( 10.D^P[7] * 2.8 * double(!csi.mH) / P[4] ) * P[5]
	Inu = ( bbc * (1 - exp(-tauc)) ) + ( bbw * (1 - exp(-tauw)) )
        return, Inu
  end
;---------------------------------------------------

;============== BEGIN GREYBODY.PRO =================

;---------------------------------------------------


function greybody, x, fluxarr, ERRARR=errarr, XUNIT=xunit, FUNIT=funit, PAR_INFO=par_info, FWHM=fwhm, $
			 PXWIDTH=pxwidth, MODELNO=modelno, KCLRCOR=kclrcor, THRESHOLD=threshold, WEIGHT=weight, $
			 PLOT=plot, SAVEPLOT=saveplot, AUTOCALL=autocall, KCOR_FILES=kcor_files

    ;;not using common blocks here because I want this function to be usable
    ;; independently as well - autocall kwarg set if called by greybody function
    ;;Use gbpstruct.pro to make par_info for direct command-line usage
    ;;AUTOCALL=1 when fxn called by greybody3D

	Nwav = N_elements(x)
	IF (Nwav LE 0) THEN RETALL, 'Error: missing wavelength array'
	IF (Nwav NE size(reform(fluxarr),/dimensions)) THEN BEGIN
	    MESSAGE,"Error: single pixel stack extraction failed" + string(7b),/INFO
	    RETALL,Nwav
          ENDIF

	IF ~keyword_set(autocall) THEN BEGIN
	    show=!null

    	    IF n_elements(xunit) LT 1 THEN xunit='um'
	    ;nudat = x2nu(x,xunit)
    	    IF (n_elements(funit) LT 1) THEN funit='Jy/beam'
            IF (strpos(funit,'beam') NE -1) AND (N_elements(fwhm) EQ 0) THEN BEGIN
               fwhm = 37.d
               print,"Data in Jy/beam but beam FWHM not supplied -- assuming Mopra resolution"
            ENDIF ELSE beam = !NULL
            IF (strpos(funit,'pixel') NE -1) AND (N_elements(pxwidth) EQ 0) THEN BEGIN
               pxw = 12.d
               print,"Data in Jy/pixel but pixel size not supplied -- assuming Mopra resolution"
            ENDIF ELSE pxw = !NULL
    	    IF (n_elements(threshold) LT 1) THEN threshold = 1e-03

	    ;fdat = fucon(fluxarr,funit,xunit,0)
            SIdata = nuFnu2SI(x,fluxarr,xunit,funit,beam=fwhm,pxsz=pxw)
            nudat = SIdata.nu
	    print, nudat
            fdat = SIdata.flux

	    IF (n_elements(errarr) NE 0) THEN BEGIN
		edat = nuFnu2SI(x,errarr,xunit,funit,beam=fwhm,pxsz=pxw,/xoff)
		weight = !null
	    ENDIF ELSE BEGIN
	    	IF (n_elements(weight) NE 0) THEN wdat = weight
	    ENDELSE
		
	    bounds = (threshold NE 0) ? nuFnu2SI(x,threshold,xunit,funit,beam=fwhm,pxsz=pxw,/xoff) : 0
	    ;; for some reason the logic doesn't work unless the bounds is exactly 0
	    ;; assumes threshold & flux are in the same units

	ENDIF ELSE BEGIN
	    ;;fwhm not needed
	    nudat = x
	    fdat = reform(fluxarr)
	    IF (n_elements(errarr) NE 0) THEN BEGIN
		edat = reform(errarr)
		weight = !null
	    ENDIF ELSE BEGIN
	    	IF (n_elements(weight) NE 0) THEN wdat = reform(weight)
	    ENDELSE
	    bounds = (n_elements(threshold) GT 1) ? reform(threshold) : threshold
	    

	ENDELSE

	IF (n_elements(errarr) EQ 0) AND (n_elements(weight) EQ 0) THEN wdat=double((fdat/max(fdat))*(max(nudat)/nudat))

	IF (n_elements(modelno) EQ 0) THEN modelno = 0
	CASE modelno OF
	    -1:model = 'test2pt' ;;-->diagnostic only, no call to MPFIT
	    0:model = 'mbb1opthin'
	    1:model = 'mbb1opthick'
	    2:model = 'mbb2comp' ;;--> probably not appropriate for Far-IR data
	    3:model = 'mbb3comp' ;;--> this one's not working yet, probably insufficient data
	    4:model = 'mbb2src' ;;--> 2 components in emission
	    ELSE:MESSAGE,"Error: invalid model number" + string(7b),/INFO
	  ENDCASE

	print,"fitting model "+model
	;;NOTE: mpfitfun does not like masked data - take care of it topside

	cc=0
	WHILE (cc LT 2) DO BEGIN
	    IF (edat NE !null) THEN BEGIN
	        parfit = mpfitfun(model, nudat, fdat, edat, parinfo=par_info, PERROR=perror, $
			COVAR=covar, STATUS=status, ERRMSG=errmsg)
	    ENDIF ELSE parfit = mpfitfun(model, nudat, fdat, parinfo=par_info, WEIGHT=wdat, $
		  		PERROR=perror, COVAR=covar, STATUS=status, ERRMSG=errmsg)
	    
	    check1 = reform(perror[where(par_info[*].fixed EQ 0)] GT 0.5*parfit[where(par_info[*].fixed EQ 0)])
	    check0 = reform(perror[where(par_info[*].fixed EQ 0)] EQ 0)
  	    IF (modelno EQ 4) THEN BEGIN
	       IF ( (status LE 0) OR (total(check0) GT 0) OR (total(check1) GT 0) OR (perror[3] GT 1.d) ) THEN BEGIN
		  print, 'Two-temp component SED fit failed - attempting 1-component fit'
		  parinfo2 = par_info[0:5]
		  ;parinfo2[0].limits[1]=[150.] ;;this may actually make things worse
		  parfit = make_array(N_elements(par_info), value=!values.d_nan)
		  perror2=perror[0:5]
		  parfit[0:5] = mpfitfun('mbb1opthick', nudat, fdat, edat, parinfo=parinfo2, PERROR=perror2, $
			COVAR=covar, STATUS=status, ERRMSG=errmsg)
	          perror[0:5]=perror2
		  perror[6:-1]=!values.d_nan
               ENDIF
            ENDIF

	    check0 = reform(perror[where(par_info[0:5].fixed EQ 0)] EQ 0)
	    checkN = finite(reform(perror[where(par_info[0:5].fixed EQ 0)]))
            IF (status LE 0) OR (total(check0) GT 0) OR (total(checkN) LT N_elements(perror[where(par_info[0:5].fixed EQ 0)])) THEN BEGIN
	    ;;only give up if single-T fit fails too
               PRINT,'status: '+string(status)+' - '+errmsg
               funcvals = make_array(N_elements(nudat), value=!values.d_nan)
	       parfit = make_array(N_elements(par_info), value=!values.d_nan)
               perror = make_array(N_elements(par_info), value=!values.d_nan)                 
	    ENDIF ELSE funcvals = call_function(model,nudat, parfit)

	    cc=cc+1
	    IF (cc GT 1) or ~keyword_set(kclrcor) THEN BREAK $
	    ELSE IF (edat NE !null) AND (total(check0) EQ 0) THEN $
		kcorrect,nudat,fdat,edat,parfit,perror,modelno=modelno,ccfiles=kcor_files


	ENDWHILE

	IF (N_elements(edat) EQ N_elements(nudat)) THEN BEGIN
	   result = {nu:nudat, flux:fdat, errflux:edat, params:parfit, perrs:perror, fmodel:funcvals}
	ENDIF ELSE result = {nu:nudat, flux:fdat, params:parfit, perrs:perror, fmodel:funcvals}
	IF keyword_set(plot) THEN sedplot,result,modelno=modelno,/listp
return, result
end
	
