function gbpstruct, T, NU0=nu0, BETA=beta, LOGN=logN, G2D=g2d, KAPPA0=kappa0, $
		TBG=Tbg, LOGNBG=logNbg, VARY_ALL=vary_all
		;;don't want to do taglist b/c it's a pain in the arse to parse
		;;taglist, PAR_INIT=parmInit, PAR_MIN=Pmin, PAR_MAX=Pmax, PAR_FIX=pfix
; PURPOSE:
;	make structure for values taken by greybody.pro & greybody_fit.pro
; CALLING:
;	gbp = gb_struct(...)
;	each arg/kwarg can be a float, or 2-, or 3-tuple:
;	arg = init_value -OR-
;	arg = [init_value, bool_fixed] -OR-
;	arg = [init_value, min, max]
;	*use -1 for min/max to switch to default min/max respectively
; INPUTS:
;	T = initial guess of temperature for greybody fit
; KEYWORDS: each can be a single value or a tuple of 2, 3, or 4 entries
;	NU0 = fiducial wavelength to fit kappa & beta w.r.t.
;	BETA = initial guess for emissivity
;	KAPPA0 = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg*
;	lOGNCOLD = initial guess for column density from cold dust, in Hmol/m^2
;	G2D = initial guess for gas-to-dust ratio; default is 133
;	TBG = initial guess for temp of 2nd component
;	lOGNBG= initial guess for column density for optional 2nd component, Hmol/m^-2
;	VARY_ALL = Boolean: if set, unfix all variables (not recommended)
; OUTPUTS:
;	PAR_INFO structure to use in MPFITFUN 
; COMMON BLOCKS:
;	TBD
; EXTERNAL CALLS:
;	greybody.pro
;	
; HISTORY:
;	Written: Rebecca Pitts, 2016.
;---------------------------------------------------
	;; Defaults
	keylist = ['T', 'NU0', 'BETA', 'LOGN', 'G2D', 'KAPPA0', 'TBG', 'LOGNBG']
	IF (n_elements(Tbg) GT 0) OR (n_elements(logNbg) GT 0) THEN BEGIN
	    nkeys = 8
	    defvals = [T, (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.55D, 130.D, 24.D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 20.D, 10.D, 0.001D, 50.D, 18.D]
	    defmaxa = [50.D, double(3.0e+13), 3.D, 32.D, 2000.D, 0.65D, 1500.D, 28.D]
	  ENDIF ELSE BEGIN
	    nkeys = 6
	    defvals = [T, (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.55D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 20.D, 10.D, 0.001D]
	    defmaxa = [50.D, double(3.0e+13), 3.D, 32.D, 2000.D, 0.65D]
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
	    IF ((N_elements(scope_varfetch(element)) GT 0) && (element NE 'T')) THEN BEGIN
		var = scope_varfetch(element)
		;;^cannot do this outside this if-statement or scope_varfetch GE 1 will always be true
		par_info[index].value = var[0]
		nel = N_elements(var)
		CASE 1 OF ;;i.e. if statements before colon are true
		    (nel EQ 1):	BREAK
		    (nel EQ 2): IF (var[1] EQ 0) THEN par_info[index].fixed = 0 ELSE var[1] = 1
		    (nel EQ 3): BEGIN
		        par_info[index].fixed = 0
		        IF (var[1] NE -1) THEN par_info[index].limits[0] = var[1]
		        IF (var[2] NE -1) THEN par_info[index].limits[1] = var[2]
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

  return, par_info
end
