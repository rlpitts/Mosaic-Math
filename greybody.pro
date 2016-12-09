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
;		in either Jy / arcsec^2, MJy/sr, or W/m^2/hz
; KEYWORDS:
;	XUNIT = string units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = alphanumeric string units of fdat
;		accepts 'MJysr', 'Jyarcsec2', or 'Wm2hz' 
;		default is 'Jyarcsec2'
;	THRESHOLD = lower limit of valid flux values
;		default is 1e-05, assuming FUNIT='Jyarcsec2'
;		will incorporate upper limits later
;	VARY_ALL = Boolean: if set, unfix all variables (not recommended)
;  ***Each of the following can be a single value or a tuple of 2 or 3 entries
;	Tcold = initial guess of temperature for greybody fit
;	NU0 = fiducial wavelength to fit kappa & beta w.r.t.
;	BETA = initial guess for emissivity
;	KAPPA0 = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg
;		(note: 1 m^2/kg = 10 cm^2/g - please do your own conversion)
;	lOGNCOLD = initial guess for column density from cold dust, in Hmol/m^2
;	G2D = initial guess for gas-to-dust ratio; default is 133
;	THOT = initial guess for temp of optional 2nd component
;	lOGNHOT = initial guess for column density for optional 2nd component, Hmol/m^-2
; OUTPUTS:
;	
; COMMON BLOCKS:
;	None b/c I'm a Python native & these things bug the scheisse out of me
;	(I kid - it's a neat trick but I'm writing most of these
;	modules to work alone OR together)
; EXTERNAL CALLS:
;	funit_conv.pro
;	gbpstruct.pro (trying to decommish that one)
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


function mbbopthin, X, P ;deriv, dP 
	;; P[0] = T, P[1] = nu_0, P[2] = beta, P[3] = log(N_col), P[4] = G2D, P[5] = kappa_0
	;; mu = mmw per unit H - assumes 71% H, 27% He, 2% Z
	;;from Arthur Cox 2000 reproduction of "Allen's Astrophysical Quantities" (1955) - everyone uses this
	;;alternative cited by Wikipedia doesn't seem to be reliable
	;; assume mu = 2.8 & omega accounted for by xunit coversion
	Inu = planck(X, P[0]) * ((X / P[1])^P[2] ) * ( 10^P[3] * double(2.8 * !const.mH) / P[4] ) * P[5]

;;	if N_params() ge 3 then begin
;;	    Pderiv = dblarr(N_elements(X),N_elements(P))
;;
;;	    dT = Inu / (P[0]^2 * (exp( double(!const.h) * X / ( double(!const.k) * P[0] ) ) - 1) )
;;	    dnu0 = -P[2] * Inu / P[1]
;;	    dbeta = alog(P[2]) * Inu
;;	    dNcol = ( alog(10)^2 ) * alog10(P[3]) * Inu
;;	    dG2D = -Inu / P[4]
;;	    dkappa = Inu/P[5]
;;
;;	    Pderiv[*,0] = dT
;;	    Pderiv[*,1] = dnu0
;;	    Pderiv[*,2] = dbeta
;;	    Pderiv[*,3] = dNcol
;;	    Pderiv[*,4] = dG2D
;;	    Pderiv[*,5] = dkappa
;;
;;	    dP = Pderiv
;;	  endif

	return, Inu
  end
;---------------------------------------------------


;; use APEX SAM threshold to determine fg/bg contamination (maybe)
;;function rt1d_mbb, X, P,
;;	;; P[0] = Tc, P[1] = Th, P[2] = nu_0, P[3] = beta, P[4] = log(N_col_c),
;;	;; P[5] = log(N_col_h), P[6] = G2D, P[7] = kappa_0
;;	;; assume mu = 2.8 & no need for omega
;;    	bb1 = planck(X, P[0])
;;	bb2 = planck(X, P[1])
;;	plaw = ((X / P[2])^P[3] )
;;	Sig1 = 10^P[4] * double(2.8 * !const.mH) / P[6]
;;	Sig2 = 10^P[5] * double(2.8 * !const.mH) / P[6]
;;	return, [ ( bb1 * (1 - exp( -plaw * Sig1 ) ) + ( bb2 * plaw * Sig2 ) * P[7] ) ]
;;  end
;---------------------------------------------------


;;function PAHmodel,...
;---------------------------------------------------


function x2nu, x, u
	CASE u OF
	    'm': nudat = reverse((double(!const.c))/x)
	    'mm': nudat = reverse((double(!const.c) * 1e+3) / x)
	    'um': nudat = reverse((double(!const.c) * 1e+6) / x)
	    'nm': nudat = reverse((double(!const.c) * 1e+9) / x)
	    'Hz': nudat = double(x)
	    'MHz': nudat = double(x) * 1e+6
	    'GHz': nudat = double(x) * 1e+9
	    'THz': nudat = double(x) * 1e+12
	    ELSE: message, 'Not a valid wavelength or frequency unit. Check spelling and case.'
          ENDCASE
  	return, nudat
  end
;---------------------------------------------------


function fucon, f, u, xu, invrs
    ;; doesn't matter what invrs is if not 0
    IF ((n_elements(f) GT 1) && (where(strmatch(xu,'*Hz',/fold_case)) EQ -1)) THEN $
     	  rf = reverse(reform(double(f))) ELSE rf = reform(double(f))
    ;; only reverse f if x was in wavelength units; above checks to see if Hz in xunit
    ;; if invrs=1, flux array should be already reversed so reversing again should set it right
    IF ((invrs eq 0) || ~keyword_set(invrs)) THEN BEGIN
	CASE u OF
	    ;; 1 MJy/sr = 2.350443e-5 Jy/arcsec^2 = 10^-20 W/m^2/Hz
	    'Jyarcsec2': newf = (rf / double(2.350443e-5)) * 1e-20
	    'MJysr' : newf = rf * 1e-20
	    'Wm2Hz' : newf = rf
	    ELSE: message, "Not a valid flux density unit. Enter 'MJysr', 'Jyarcsec2', or 'Wm2Hz'."
	  ENDCASE
      ENDIF ELSE BEGIN
	CASE u OF
	    'Jyarcsec2': newf = (rf * double(2.350443e-5)) * 1e+20
	    'MJysr' : newf = rf * 1e+20
	    'Wm2Hz' : newf = rf
	    ELSE: message, "Not a valid flux density unit. Enter 'MJysr', 'Jyarcsec2', or 'Wm2Hz'."
	  ENDCASE
      ENDELSE
  return, newf
  end
;---------------------------------------------------


function greybody, x, fluxarr, XUNIT=xunit, FUNIT=funit, PAR_INFO=par_info, $
			THRESHOLD=threshold, WEIGHT=weight, AUTOCALL=autocall ;;,$
			;;PAR_INIT=parmInit, PAR_MIN=Pmin, PAR_MAX=Pmax, PAR_FIT=parmFit, $
			;;NPARFIX=nparfix, PLOTFIT=plotfit, SHOW=show, MAXITER=maxit, $
			;;INFO=info, VERBOSE=verbose

    ;;not using common blocks here because I want this function to be usable
    ;; independently as well - autocall kwarg set if called by greybody function
    ;;Use gbpstruct.pro to make par_info for direct command-line usage

	Nwav = N_elements(x)
	IF (Nwav LE 0) THEN RETALL, 'Error: missing wavelength array'
	IF (Nwav NE size(reform(fluxarr),/dimensions)) THEN BEGIN
	    MESSAGE,"Error: single pixel stack extraction failed" + string(7b),/INFO
	    RETALL,Nwav
          ENDIF

	IF ~keyword_set(autocall) THEN BEGIN
    	    IF n_elements(xunit) LT 1 THEN xunit='um'
	    nudat = x2nu(x,xunit)
    	    IF (n_elements(funit) LT 1) THEN funit='Jyarcsec2'
    	    IF (n_elements(threshold) LT 1) THEN threshold = 1e-05
	    fdat = fucon(fluxarr,funit,xunit,0)
	    floor = (threshold NE 0) ? fucon(threshold,funit,xunit,0) : 0
	    ;;assumes threshold & flux are in the same units
	  ENDIF ELSE BEGIN
	    nudat = x
	    fdat = reform(fluxarr)
	    floor = threshold
	  ENDELSE

	IF (n_elements(weight) EQ 0) then weight=double((fdat/max(fdat))*(max(nudat)/nudat))

	;;NOTE: mpfitfun does not like masked data - take care of it topside
	IF (n_elements(par_info) EQ 6) THEN BEGIN
	    ;;parfit = par_info[*].value[*]
	    parfit = mpfitfun('mbbopthin', nudat, fdat, parinfo=par_info, weights=weight,$
		STATUS=status, ERRMSG=errmsg)
  	    IF status LE 0 THEN message, status, errmsg
	    funcvals = mbbopthin(nudat, parfit)
	    ;print, funcvals
	  ENDIF

	IF ~keyword_set(autocall) THEN result = {nu:nudat, fdata:fdat, params:parfit, fmodel:funcvals} ELSE $
	result = parfit

return, result
end
	
