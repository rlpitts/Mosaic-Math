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
;	VARY_ALL = Boolean: if set, unfix all variables (not recommended)
;	AUTOCALL = Boolean switch to distinguish solo calls from calls by greybody3D.pro
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
;	*Assumes !const is defined, is in SI units, & contains mass of hydrogen atom
; HISTORY:
;	Written: Rebecca Pitts, 2016.
;---------------------------------------------------

function planck, nu, T
    	hc2 =  2 * double(!const.h) / double(!const.c)^2
    	hkt = double(!const.h) / ( double(!const.k) * T )
    	bb = hc2 * nu^3 / (exp( hkt * nu ) - 1)
	return, bb
  end
;---------------------------------------------------


function mbb1opthin, X, P
	;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, P[5] = kappa_0
	;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
	;;from Arthur Cox 2000 reproduction of "Allen's Astrophysical Quantities" (1955) - everyone uses this
	;;alternative cited by Wikipedia doesn't seem to be reliable
	;; assume mu = 2.8 & omega accounted for by xunit coversion

	Inu = planck(X, P[0]) * ( (X / P[1])^P[2] ) * ( 10.D^P[3] * double(2.8 * !const.mH) / P[4] ) * P[5]
	return, Inu
  end
;---------------------------------------------------


;; ASSUMPTIONS: negligible background - all emission local to source & self-attenuating
function mbb1opthick, X, P
	;; P[0] = Tc, P[1] = nu_0, P[2] = beta, P[3] = log(N_col_c), P[4] = G2D, P[5] = kappa_0
	;; assume mu = 2.8 & no need for omega
	;; PARAMETERS CANNOT BE REARRANGED (yes, it sucks)

	bb = planck(X, P[0])
	tau = ( (X / P[1])^P[2] ) * ( 10.D^P[3] * double(2.8 * !const.mH) / P[4] ) * P[5]
    	Inu = bb * ( 1 - exp(-tau) )
	return, Inu
  end
;---------------------------------------------------


function mbb2comp, X, P
        ;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, 
        ;; P[5] = kappa_0, P[6] = T_bg
        ;; (see e.g. Magnum & Shirley 2015 review)
        ;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
        ;; assume mu = 2.8, UNIFORM DUST COMPOSITION* & omega accounted for by xunit coversion
        ;;*This is a whopper if there ever was one

	bb = planck(X, P[0])
	tau = ((X / P[1])^P[2] ) * ( 10.D^P[3] * double(2.8 * !const.mH) / P[4] ) * P[5]
	bb_bg = planck(X, P[6])
	Inu = ( bb * (1 - exp(-tau)) ) + ( bb_bg * exp(-tau) )
        return, Inu
  end
;---------------------------------------------------


function mbb3comp, X, P ;deriv, dP 
;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, 
;; P[5] = kappa_0, P[6] = T_bg, P[7] = T_bg
;; (see e.g. Magnum & Shirley 2015 review)
;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
;; assume mu = 2.8, UNIFORM DUST COMPOSITION* & omega accounted for by xunit coversion
;;*again, majorly suspect assumption

	bb = planck(X, P[0])
	tau = ((X / P[1])^P[2] ) * ( 10.^P[3] * double(2.8 * !const.mH) / P[4] ) * P[5]
	bb_bg = planck(X, P[6])
	bb_fg = planck(X, P[7])
	Inu = ( bb * (1 - exp(-tau)) ) + ( bb_bg * exp(-tau) ) + bb_fg
	return, Inu
  end
;---------------------------------------------------

;============== BEGIN GREYBODY.PRO =================

;---------------------------------------------------


function greybody, x, fluxarr, ERRARR=errarr, XUNIT=xunit, FUNIT=funit, PAR_INFO=par_info,$
			 FWHM=fwhm, PXWIDTH=pxwidth, MODELNO=modelno, THRESHOLD=threshold, WEIGHT=weight,$
			 SAVEPLOT=saveplot, AUTOCALL=autocall

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

	IF (n_elements(modelno) EQ 0) THEN modelno = 1
	CASE modelno OF
	    1:model = 'mbb1opthin'
	    2:model = 'mbb1opthick'
	    3:model = 'mbb2comp' ;;--> kinda works in a different form (see iter_gb3d)
	    4:model = 'mbb3comp' ;;--> this one's not working yet, probably insufficient data
	    ELSE:MESSAGE,"Error: invalid model number" + string(7b),/INFO
	  ENDCASE

	;;NOTE: mpfitfun does not like masked data - take care of it topside
	IF (modelno LT 3) THEN BEGIN
	    IF (edat NE !null) THEN BEGIN
	        parfit = mpfitfun(model, nudat, fdat, edat, parinfo=par_info, PERROR=perror, $
			COVAR=covar, STATUS=status, ERRMSG=errmsg)
	    ENDIF ELSE parfit = mpfitfun(model, nudat, fdat, parinfo=par_info, WEIGHT=wdat, $
		  		PERROR=perror, COVAR=covar, STATUS=status, ERRMSG=errmsg)
  	    IF (status LE 0) THEN BEGIN
		MESSAGE,'status: '+str(status)+' - '+errmsg
		funcvals = make_array(N_elements(nudat), value=!values.d_nan)
		perror = make_array(N_elements(par_info), value=!values.d_nan)
	    ENDIF ELSE funcvals = mbb1opthin(nudat, parfit)
	    ;;
	  ENDIF

	IF ~keyword_set(autocall) THEN BEGIN
	    IF (N_elements(edat) EQ N_elements(nudat)) THEN BEGIN
		result = {nu:nudat, fdata:fdat, ferr:edat, params:parfit, parerrs:perror, fmodel:funcvals}
	    ENDIF ELSE result = {nu:nudat, fdata:fdat, params:parfit, parerrs:perror, fmodel:funcvals}
	ENDIF ELSE BEGIN
	    result = {params:parfit, parerrs:perror}
	ENDELSE
	;;call sedplot.pro if you want to plot
return, result
end
	
