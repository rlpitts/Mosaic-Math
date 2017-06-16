function gbpstruct, TC, NU0=nu0, BETA=beta, LOGNC=logNc, G2D=g2d, KAPPA0=kappa0, $
		TBG=Tbg, TFG=Tfg, VARY_ALL=vary_all
		;;don't want to do taglist b/c it's a pain in the arse to parse
		;;taglist, PAR_INIT=parmInit, PAR_MIN=Pmin, PAR_MAX=Pmax, PAR_FIX=pfix
; PURPOSE:
;	make structure for values taken by greybody.pro & greybody_fit.pro
; CALLING:
;	gbp = gb_struct(...)
;	each arg/kwarg can be a float, or 2-, or 3-tuple:
;	arg = init_value* -OR-
;	arg = [init_value, bool_fixed]* -OR-
;	arg = [init_value, min, max]**
;	*arg = init_value is equivalent to arg = [init_value, 1];
;	     use arg = [init_value, 0] to vary w/in default limits
;	**use -1 for min/max to switch to default min/max respectively
; INPUTS:
;	TC = initial guess of temperature for greybody fit
; KEYWORDS: each can be a single value or a tuple of 2 or 3 entries
;	NU0 = fiducial wavelength to fit kappa & beta w.r.t.
;	BETA* = initial guess for emissivity; default is 1.8
;		Other recommended values: 1.5, 2.0, 2.1, & 2.2
;		*WARNING: DEGENERATE WITH TEMPERATURE; also, must be
;		consistent with kappa0 (check e.g. Bianchi 2013)
;	KAPPA0** = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg*
;		*(note: 1 m^2/kg = 10 cm^2/g - please convert beforehand)
;		**depends on G2D: tau/N = kappa*m_H/G2D
;		At 250 um:
;		0.68 = (tau/N)*(G2D/m_H) = 9.2e-26*(124./m_H) (beta=1.8, solar neighborhood)
;		0.89 = (tau/N)*(G2D/m_H) = 9.2e-26*(162./m_H) (beta=1.8, solar neighborhood)
;		0.47, 0.51 = (tau/N)*(G2D/m_H) = (4.9, 5.3)e-26*(162./m_H) (HeViCS HLGC,beta=1.53, 1.6)
;		0.36, 0.39 = (tau/N)*(G2D/m_H) = (4.9, 5.3)e-26*(124./m_H) (HeViCS HLGC,beta=1.53, 1.6)
;	LOGNC = initial guess for column density from cold dust, in Hmol/m^2
;	G2D* = initial guess for gas-to-dust ratio; default is 124 (MW).
;		from Draine and Li 2001; other recommended values (for MW):
;		133 (Compiegne+2010 - caution: inconsistent opacity & beta)
;		162 (average from Zubko, Dwek, & Arendt 2004)
;		*WARNING: DEGENERATE WITH COLUMN DENSITY
;	TBG = initial guess for temp of 2nd component
;	TFG = initial guess for temp of 3rd component
;	VARY_ALL = Boolean: if set, unfix all variables (not recommended)
; OUTPUTS:
;	PAR_INFO structure to use in MPFITFUN 
; COMMON BLOCKS:
;	keys - contains keylist
; EXTERNAL CALLS:
;	greybody.pro
;	greybody3D.pro
; HISTORY:
;	Written: Rebecca Pitts, 2016.
;---------------------------------------------------

    common keys, keylist, par_info

	;; Defaults
	keylist = ['TC', 'NU0', 'BETA', 'LOGNC', 'G2D', 'KAPPA0', 'TBG', 'TFG']
	IF (n_elements(Tbg) GT 0) OR (n_elements(logNbg) GT 0) THEN BEGIN
	    nkeys = 8
	    defvals = [double(Tc), (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.68D, 130.D, 90.D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 22.D, 5.D, 0.001D, 50.D, 20.D]
	    defmaxa = [90.D, double(3.0e+13), 3.D, 28.D, 2000.D, 1.5D, 1500.D, 1500.D]
	  ENDIF ELSE BEGIN
	    nkeys = 6
	    defvals = [double(Tc), (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.68D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 22.D, 5.D, 0.001D]
	    defmaxa = [90.D, double(3.0e+13), 3.D, 28.D, 2000.D, 1.5D]
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
	    IF ((N_elements(scope_varfetch(element)) GT 0) && (element NE 'TC')) THEN BEGIN
		var = scope_varfetch(element)
		;;^cannot do this outside this if-statement or scope_varfetch GE 1 will always be true
		par_info[index].value = double(var[0])
		nel = N_elements(var)
		CASE 1 OF
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
	;; Models fit gas-to-dust ratio, not dust-to-gas
	IF (par_info[4].value LT 1.D) THEN par_info[4].value = 1/par_info[4].value
	IF (par_info[4].limits[0] LT 1.D) THEN par_info[4].limits[0] = 1/par_info[4].limits[0]
	IF (par_info[4].limits[1] LT 1.D) THEN par_info[4].limits[1] = 1/par_info[4].limits[1]

	;; I know this segment needs more finessing but I still don't know all the constraints
	IF ((par_info[2].value GT 1.85) && (par_info[2].value LE 2.2)) THEN BEGIN
	    ;;Second default to Herschel values
	    ;;(Magrini et al 2011, Davies et al 2012, Smith et al 2012, Bianchi 2013)
	    CASE 1 OF ;;means if any of the following is true:
	        (n_elements(kappa0) EQ 0): BEGIN
		    par_info[5].value = 0.192
		    par_info[1].value = (double(!const.c) * 1e6)/350.D
		    ;;Second default to Herschel values
		    ;;(Magrini et al 2011, Davies et al 2012, Smith et al 2012, Bianchi 2013)
		  END
		((kappa0[0] GE 0.001) && (kappa0[0] LT 0.1)): par_info[1].value = (double(!const.c) * 1e6)/850.D
		((kappa0[0] GE 0.1) && (kappa0[0] LT 0.3)): par_info[1].value = (double(!const.c) * 1e6)/350.D
		((kappa0[0] GE 0.3) && (kappa0[0] LT 0.8)): par_info[1].value = (double(!const.c) * 1e6)/250.D
		ELSE: STOP, 'opacity incompatible with built-in emissivity and fiducial frequency values.'
	      ENDCASE	
	  ENDIF

  return, par_info
end
