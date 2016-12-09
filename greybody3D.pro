;+
; NAME:
;	greybody3D.pro
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
;	greybody.pro
;	
;	gbpstruct.pro (trying to decommish that one)
;	*Assumes !const is defined, is in SI units, & contains mass of hydrogen atom
; HISTORY:
;	Written: Rebecca Pitts, 2016.
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
    ndims = size(f,/n_dimensions)
    IF ((ndims GT 0) && (where(strmatch(xu,'*Hz',/fold_case)) EQ -1)) THEN $
     	  rf = reverse(reform(double(f)),ndims) ELSE rf = reform(double(f))
    ;;reverse(a,ndims) reverses array depthwise/along its last dimension
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

function greybody3D, x, imcube, XUNIT=xunit, FUNIT=funit, TCOLD=Tcold, NU0=nu0, $
		BETA=beta, LOGNCOLD=logNcold, G2D=g2d, KAPPA0=kappa0, $
		THOT=Thot, LOGNHOT=logNhot, VARY_ALL=vary_all, THRESHOLD=threshold, $
		WEIGHT=weight, SAVEP=savep;;VERBOSE=verbose, INFO=info, SHOW=show, MAXITER=maxit

    ;;common USER_INPUT x, xunit, funit, par_info
    ;;common IR_spectrum_fit, fitControl, WghtControl
    ;;common IR_spectrum_fit0, fitInitParams
    ;;common IR_spectrum_fit1, fitPlotOptions
    ;;common IR_spectrum_fit2, FitParLimits

	Nwav = N_elements(x)
	;;/dimensions gives size along each axis, /n_dimensions counts total dimensions
	imd = size(imcube,/dimensions) ;;imd[-1] will be the same whether 1D or 3D - I checked
	IF (Nwav LE 0) THEN RETALL, 'Error: Missing wavelength array'

	IF (Nwav NE imd[-1]) THEN BEGIN
	    MESSAGE,"Error: unequal image cube depth and length of wavelength array" + string(7b),/INFO
	    RETALL,Nwav,imd
          ENDIF

    	if n_elements(xunit) lt 1 then xunit='um'
    	if n_elements(funit) lt 1 then funit='Jyarcsec2'
	nudat = x2nu(x,xunit) ;;convert flux in greybody_fit
	fcube = fucon(imcube,funit,xunit,0)
	IF n_elements(threshold) EQ 0 THEN threshold = 0.D
	floor = (threshold NE 0) ? fucon(threshold,funit,xunit,0) : 0.D
	;;mask cube before fucon - afterward, the precision needed to evaluate
	;; relational operators is too small

	;; Defaults
	keylist = ['TCOLD', 'NU0', 'BETA', 'LOGNCOLD', 'G2D', 'KAPPA0', 'THOT', 'LOGNHOT']
	IF (n_elements(THOT) GT 0) OR (n_elements(LOGNHOT) GT 0) THEN BEGIN
	    nkeys = 8
	    defvals = [15.D, (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.55D, 130.D, 24.D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 20.D, 10.D, 0.001D, 50.D, 18.D]
	    defmaxa = [50.D, double(3.0e+13), 3.D, 32.D, 2000.D, 0.65D, 1500.D, 28.D]
	  ENDIF ELSE BEGIN
	    nkeys = 6
	    defvals = [20.D, (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.55D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 20.D, 10.D, 0.001D]
	    defmaxa = [90.D, double(3.0e+13), 3.D, 32.D, 2000.D, 0.65D]
	  ENDELSE

	par_info = replicate({value:0.D, fixed:1, limited:[1,1], limits:[0.D,0]}, nkeys)
	par_info[*].value = defvals
	par_info[*].limits[0] = defmina
	par_info[*].limits[1] = defmaxa
	par_info[0].fixed = 0
	
	IF keyword_set(vary_all) EQ 1 THEN par_info[*].fixed = 0
	
	;; Command-line overrides
	;; NOTE: scope_varfetch() looks for kwarg names to the right of the '=', not the left O_o
	FOREACH element, keylist, index DO BEGIN
	    IF (N_elements(scope_varfetch(element)) GT 0) THEN BEGIN
		var = scope_varfetch(element)
		;;^can't do outside this if-statement or (scope_varfetch GE 1) will always be true
		par_info[index].value = double(var[0])
		nel = N_elements(var)
		CASE 1 OF ;;i.e. if statements before colon are true
		    (nel EQ 1):	BREAK
		    (nel EQ 2): IF (var[1] EQ 0) THEN par_info[index].fixed = 0 ELSE var[1] = 1
		    (nel EQ 3): BEGIN
		        par_info[index].fixed = 0
		        IF (var[1] NE -1) THEN par_info[index].limits[0] = double(var[1])
		        IF (var[2] NE -1) THEN par_info[index].limits[1] = double(var[2])
		      END
		    (nel GT 3): STOP, 'Each kwarg takes a 3-tuple at most. Check '+element
		  ENDCASE
	      ENDIF
	  ENDFOREACH

	;; Corrective overrides of the overrides (Big Brother is watching)
	IF (par_info[4].value LT 1.D) THEN par_info[4].value = 1/par_info[4].value
	IF (par_info[4].limits[0] LT 1.D) THEN par_info[4].limits[0] = 1/par_info[4].limits[0]
	IF (par_info[4].limits[1] LT 1.D) THEN par_info[4].limits[1] = 1/par_info[4].limits[1]

	;; I know this segment needs more finessing but I still don't know all the constraints
	IF ((par_info[2].value GT 1.9) && (par_info[2].value LT 2.1)) THEN BEGIN
	    ;;Second default to Herschel values
	    ;;(Magrini et al 2011, Davies et al 2012, Smith et al 2012)
	    CASE 1 OF ;;means if any of the following is true:
	        (n_elements(kappa0) EQ 0): BEGIN
		    par_info[5].value = 0.192
		    par_info[1].value = (double(!const.c) * 1e+6)/350
		  END
		((kappa0[0] GE 0.001) && (kappa0[0] LT 0.1)): par_info[1].value = (double(!const.c) * 1e+6)/850
		((kappa0[0] GE 0.1) && (kappa0[0] LT 0.4)): par_info[1].value = (double(!const.c) * 1e+6)/350
		((kappa0[0] GE 0.4) && (kappa0[0] LT 0.6)): par_info[1].value = (double(!const.c) * 1e+6)/250
		ELSE: STOP, 'opacity incompatible with built-in emissivity and fiducial frequency values.'
	      ENDCASE	
	  ENDIF

	;;mask bad/masked data - NOTE: values are 'good' where mask is TRUE
	;;print, size(fcube,/dimensions),fcube[20,20,3]
	mask = fcube[*,*,*] GT floor
	mfcube = fcube*mask
	goodijk = where(mask[*,*,*],count,/null)
	print, floor, count, min(fcube[goodijk])
	;;pkeys = keylist[where(par_info.fixed[*] EQ 0, /null)]
	freepars = where(par_info[*].fixed EQ 0, /null)
	parcube = fltarr(imd[0],imd[1],nkeys);;n_elements(freepars)
	;;I'll work out how to make it a hash later
	;;Our crusty old version of IDL doesn't do dictionaries or HasValues >:[

	IF N_elements(weight) EQ 0 THEN weight=!null

	FOR I=0,imd[0]-1 DO BEGIN
	    FOR J=0,imd[1]-1 DO BEGIN
		spec = reform(fcube[I,J,*])
		goodi = where(spec GT floor, /null)
		fatali = where(spec[1:-2] LT floor, count, /null)
		;print, (goodi EQ !null) ? spec : goodi
		CASE 1 OF
		    (count GT 0): parcube[I,J,*] = !values.d_nan
	    	    (N_elements(goodi) LT N_elements(par_info)): parcube[I,J,*] = !values.d_nan
		    ((spec[0] LT floor) || (spec[-1] LT floor)): BEGIN
			parcube[I,J,*] = reform(greybody(nudat[goodi], fcube[I,J,goodi], PAR_INFO=par_info,$
					   THRESHOLD=floor, WEIGHT=weight[goodi], /autocall))
		      END
		    ELSE: BEGIN
			parcube[I,J,*] = reform(greybody(nudat, fcube[I,J,*], PAR_INFO=par_info,$
					   THRESHOLD=floor, WEIGHT=weight, /autocall))
		      END
		  ENDCASE
	      ENDFOR
	  ENDFOR
	results = {nu:nudat,flux:fcube,params:parcube}
	;;IF keyword_set(verbose) THEN 
	IF keyword_set(savep) THEN save, /variables, VERBOSE=verbose, FILENAME = 'parcube.sav'

;;How to plot (don't use contour!) - e.g. for column density
;; mask=where3d(parcube.flux[*,*,-1] GT 0,xind=x,yind=y)
;; Narr=parcube.params[*,*,3]
;; Ncmap=image(Narr*mask,TITLE='Region 26 H_2 Column Density Map', RGB_TABLE=15)
;; cb = colorbar(TARGET=Ncmap,title='H_2 Column Density (m^-2)')

    return, results
end
	
