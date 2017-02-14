function greybody3D, x, imcube, ERRCUBE=errcube, XUNIT=xunit, FUNIT=funit, FWHM=fwhm,$
		TC=Tc, NU0=nu0, BETA=beta, LOGNC=logNc, G2D=g2d, KAPPA0=kappa0, $
		TBG=Tbg, LOGNBG=logNbg, VARY_ALL=vary_all, THRESHOLD=threshold, $
		WEIGHT=weight, MODELNO=modelno, SAVEP=savep, PLOTMAP=plotmap ;;, INFO=info
		;;VERBOSE=verbose, MAXITER=maxit
;+
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
;	gb = greybody3D(x, fdat, ... see kwarg list)
;	***Some kwargs can be a float, or 2-, or 3-tuple:
;	arg = init_value -OR-
;	arg = [init_value, bool_fixed] -OR-
;	arg = [init_value, min, max]
;	      (use -1 for min/max to sub in default min/max respectively)
; INPUTS:
;	x = array of wavelengths or frequencies, in any units allowed for XUNIT.
;		defaults to wavelength in microns
;	fluxarr = 3D array of fluxes, i.e. stack of images, at each wavelength
;		in either Jy/asec^2, MJy/sr, Jy/beam, or W/m^2/hz
; KEYWORDS:
;	ERRCUBE = 1D or 3D array of errors in flux, same units as fluxarr;
;		if 1D, must be same length as x
;		overrides weight if both are present
;	WEIGHT = float array of weights for each flux, nominally in 
;		units of 1/flux; nullified if errcube present
;		see greybody.pro for default
;	XUNIT = string units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = string units of fdat
;		accepts 'MJy/sr','Jy/beam','Jy/asec2', 'W/m2/Hz', or all of the
;		above without the '/'
;	THRESHOLD = lower limit of valid flux values
;		default is 1e-05, assuming FUNIT='Jyarcsec2'
;		will incorporate upper limits later
;	VARY_ALL = Boolean: if set, unfix all variables (NOT RECOMMENDED - usu.
;		insufficient data to constrain fit, plus Beta & T are degenerate)
;	SAVE_P = Boolean: if set, save cube of fitting parameters & errors to
;		parcube.sav (IDL-accessible-only .sav file)
;  ***Each of the following can be a single value or a tuple of 2 or 3 entries
;	TC = initial guess of temperature for greybody fit
;	NU0 = fiducial wavelength to fit kappa & beta w.r.t.
;	BETA = initial guess for emissivity
;	KAPPA0 = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg
;		(note: 1 m^2/kg = 10 cm^2/g - please do your own conversion)
;	LOGNC = initial guess for column density from dust at clump/core/target
;		in H2 molecules * m^-2
;	G2D = initial guess for gas-to-dust ratio; default is 133
;	TBG = initial guess for temp of optional 2nd component
;	LOGNBG = initial guess for column density for optional 2nd component, H2mol/m^-2
;		 in H2 molecules * m^-2
; OUTPUTS:
;	RESULTS: a struct containing cubes of fitting params & their errors @ each pixel
;	Also can save parameter maps & error maps - still working on 
; COMMON BLOCKS:
;	None b/c I'm a Python native & these things bug the scheisse out of me
;	(I kid - it's a neat trick but I'm writing most of these
;	modules to work alone OR together)
; EXTERNAL CALLS:
;	nuFnu2SI.pro
;	greybody.pro & its helper functions
;	*Assumes !const is defined, is in SI units, & contains mass of hydrogen atom
; HISTORY:
;	Written: Rebecca Pitts, 2016.
;	Jan 2017 - updated to run with nuFnu2SI.pro
;		helper fxns x2nu and fucon decomissioned
;	Feb 2017 - added error cube & weight options, plotmap routine
;
;---------------------------------------------------
    ;;common USER_INPUT x, xunit, funit, par_info
    ;;common IR_spectrum_fit, fitControl, WghtControl
    ;;common IR_spectrum_fit0, fitInitParams
    ;;common IR_spectrum_fit1, fitPlotOptions
    ;;common IR_spectrum_fit2, FitParLimits

	Nwav = N_elements(x)
	IF (Nwav LE 0) OR (N_elements(imcube) LE 0) THEN RETALL, 'FatalError: Missing required input (check documentation)'
	;;/dimensions gives size along each axis, /n_dimensions counts total dimensions
	imd = size(imcube,/dimensions) ;;imd[-1] will be the same whether 1D or 3D - I checked
	IF (Nwav NE imd[-1]) THEN BEGIN
	    MESSAGE,"Error: unequal image cube depth and length of wavelength array" + string(7b),/INFO
	    RETALL,Nwav,imd
         ENDIF
        IF (N_elements(fwhm) LE 0) THEN fwhm = 37.0 ;;Mopra resolution

    	IF n_elements(xunit) LT 1 THEN xunit='um'
    	IF n_elements(funit) LT 1 THEN funit='Jybeam'
        SIdata = nuFnu2SI(x,imcube,xunit,funit,beam=fwhm)
        nudat = SIdata.nu
        fcube = SIdata.flux

	IF (N_elements(errcube) GT 0) THEN BEGIN
	    erd = size(errcube,/dimensions)
	    IF (erd[-1] NE imd[-1]) THEN BEGIN
	    	MESSAGE,"Error: image cube and error cube dimensions do not match" + string(7b),/INFO
	    	RETALL,imd,erd
	    ENDIF
	    
	    ecube = nuFnu2SI(x,errcube,xunit,funit,beam=fwhm,/xoff)
	    IF (N_elements(weight) GT 0) THEN BEGIN
		PRINT, 'Error array takes precedence over weight array; weight nullified'
		weight=!null
	    ENDIF
	ENDIF ELSE IF (N_elements(errcube) EQ 0) AND (N_elements(weight) GT 0) THEN BEGIN
	    wrd = size(weight,/dimensions)
	    IF (wrd[-1] NE imd[-1]) THEN BEGIN
	    	MESSAGE,"Error: image cube and weight cube dimensions do not match" + string(7b),/INFO
	    	RETALL,imd,wrd
	    ENDIF
	    ecube = !null
	ENDIF ELSE IF (N_elements(errcube) EQ 0) AND (N_elements(weight) EQ 0) THEN BEGIN
	    weight = !null
	    ecube = !null
	    PRINT, 'Warning: defaulting to internal peak-biased weight scheme.'
	    PRINT, 'Do not trust chi-squared values or parameter errors, if they print.'
	ENDIF

	IF n_elements(threshold) EQ 0 THEN threshold = 0.D
	bounds = (threshold NE 0) ? nuFnu2SI(x,threshold,xunit,funit,beam=fwhm,/xoff) : 0.D
	;; has precision issues when bounds =/= 0

	;; Defaults
	keylist = ['TC', 'NU0', 'BETA', 'LOGNC', 'G2D', 'KAPPA0', 'TBG', 'LOGNBG']
	IF (n_elements(TBG) GT 0) OR (n_elements(LOGNBG) GT 0) THEN BEGIN
	    nkeys = 8
	    modelno = 2 ;;
	    defvals = [15.D, (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.55D, 130.D, 24.D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 20.D, 10.D, 0.001D, 50.D, 18.D]
	    defmaxa = [50.D, double(3.0e+13), 3.D, 32.D, 2000.D, 0.65D, 1500.D, 28.D]
	  ENDIF ELSE BEGIN
	    nkeys = 6
	    modelno = 1
	    defvals = [20.D, (double(!const.c) * 1e+6)/250.D, 1.8, 26.D, 124.D, 0.55D]
	    defmina = [2.73D, double(3.0e+11), 0.9D, 20.D, 10.D, 0.001D]
	    defmaxa = [90.D, double(3.0e+13), 3.D, 32.D, 2000.D, 0.65D]
	  ENDELSE

	;;IF (n_elements(modelno) EQ 0) THEN modelno = 1
	CASE modelno OF
	    1:model = 'mbb1opthin'
	    2:model = 'mbb2opthin'
	    ;;3:'mbbfrt' -- this one's not ready yet
	    ELSE:MESSAGE,"Error: invalid model number" + string(7b),/INFO
	ENDCASE

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
	mask = fcube[*,*,*] GT bounds
	mfcube = fcube*mask
	goodijk = where(mask[*,*,*],count,/null)
	;;print, bounds, count, min(fcube[goodijk])
	;;pkeys = keylist[where(par_info.fixed[*] EQ 0, /null)]
	freepars = where(par_info[*].fixed EQ 0, count, /null)
	dof=count
	parcube = fltarr(imd[0],imd[1],nkeys);;n_elements(freepars)
	puncube = fltarr(imd[0],imd[1],nkeys);;n_elements(freepars)
	;;I'll work out how to make it a hash later
	;;Our crusty old version of IDL doesn't do dictionaries or HasValues >:[

	print, 'par_info: ', par_info
	print, 'dof: ', dof

	FOR I=0,imd[0]-1 DO BEGIN
	    FOR J=0,imd[1]-1 DO BEGIN
		spec = reform(fcube[I,J,*])
		IF (ecube NE !null) THEN BEGIN
		    IF (N_elements(ecube) LE Nwav) THEN err = reform(ecube) ELSE $
			err = reform(ecube[I,J,*])
		ENDIF ELSE err = !null
		IF (weight NE !null) THEN BEGIN
		    IF (N_elements(weight) LE Nwav) THEN wt = reform(weight) ELSE $
			wt = reform(weight[I,J,*])
		ENDIF ELSE wt = !null
		IF (err NE !null) THEN goodi = where(spec GT err, /null) ELSE goodi = where(spec GT bounds, /null)
		IF (err NE !null) THEN fatali = where(spec[1:-2] LE err[1:-2]/10, count, /null) ELSE $
		    fatali = where(spec[1:-2] LE bounds, count, /null)
		IF (count GT 0) OR (N_elements(goodi) LT N_elements(par_info)) THEN BEGIN
		    parcube[I,J,*] = !values.d_nan
		    puncube[I,J,*] = !values.d_nan
	    	ENDIF ELSE IF ((spec[0] LT bounds) || (spec[-1] LT bounds)) THEN BEGIN
		    IF (N_elements(err) LE 1) THEN errarr = !null ELSE errarr = err[goodi]
		    IF (N_elements(wt) LE 1) THEN wtarr = !null ELSE wtarr = wt[goodi]
		    print, N_elements(nudat[goodi]), N_elements(fcube[I,J,goodi]), N_elements(errarr[goodi])
		    ;;I could use spec here instead of fcube[I,J,*], but I prefer not to
		    gbstruct = greybody(nudat[goodi], reform(fcube[I,J,goodi]), ERRARR=errarr, FWHM=fwhm, $
				PAR_INFO=par_info, MODELNO=modelno, THRESHOLD=bounds, WEIGHT=wtarr, /autocall)
		    parcube[I,J,*] = reform(gbstruct.params)
		    puncube[I,J,*] = reform(gbstruct.parerrs)
		ENDIF ELSE BEGIN
		    gbstruct = greybody(nudat, reform(fcube[I,J,*]), ERRARR=err, FWHM=fwhm, $
				PAR_INFO=par_info, MODELNO=modelno, THRESHOLD=bounds, WEIGHT=wt, /autocall)
		    parcube[I,J,*] = reform(gbstruct.params)
		    puncube[I,J,*] = reform(gbstruct.parerrs)
		ENDELSE
	    ENDFOR
	ENDFOR
	results = {nu:nudat,flux:fcube,params:parcube,perrs:puncube}
	;;IF keyword_set(verbose) THEN 
	IF keyword_set(savep) THEN save, /variables, VERBOSE=verbose, FILENAME = 'parcube.sav'

	IF keyword_set(plotmap) THEN BEGIN
	    Narr=results.params[*,*,3]
	    Nerr=results.perrs[*,*,3]
	    Narr[where3d(results.flux[*,*,-1] LE 0,xind=x,yind=y)]=!values.d_nan
	    Nerr[where3d(results.flux[*,*,-1] LE 0,xind=x,yind=y)]=!values.d_nan

	    Ncmap=image(Narr,TITLE='Region 26 $H_2$ Column Density Map', RGB_TABLE=15) ;;need wildcard regions
	    cb = colorbar(TARGET=Ncmap,title='log$N(H_2)$ ($m^{-2}$)', orientation=1, POSITION=[0.2,0.2,0.25,0.8])
	    ;; need user input filenames & dir names here, but I'll figure that out later
	    Ncmap.save,"~/mosaic/scratch/NH2map_R26_2dof.png", border=10, resolution=400, /close
	    print,"NH2map_R26_2dof.png saved to ~/mosaic/scratch/"

	    Ncerrmap=image(Nerr,TITLE='Region 26 Error map for fitted $H_2$ Column Density', RGB_TABLE=15) ;;need wildcard regions
	    cb = colorbar(TARGET=Ncerrmap,title='$\sigma$(log$N(H_2)$) ($m^{-2}$)', orientation=1, POSITION=[0.2,0.2,0.25,0.8])
	    ;; need user input filenames & dir names here, but I'll figure that out later
	    Ncerrmap.save,"~/mosaic/scratch/NH2map_err_R26_2dof.png", border=10, resolution=400, /close
	    print,"NH2map_err_R26_2dof.png saved to ~/mosaic/scratch/"

	ENDIF

    return, results
end
	
