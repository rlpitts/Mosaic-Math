;+
; NAME:
;	greybody.pro
; PURPOSE:
;	Compute and fit greybody emission spectrum to given flux data
;	in units of Janskys/arcsec^2 at for an array of wavelengths.
;	Uses the following relationships:
;		kappa = kappa_0*(wl_0/wl)^beta
;		kappa_0*N_H = tau_0, where N = col. density & tau << 1
;	If tau = absorption cross-section & N_H = # of H atoms, then
;		kappa*mu*m_H/G2D = tau/N_H = (tau_0/N_H)*(wl_0/wl)^beta
;	     -->tau = N_H*kappa_0*(mu*m_H/G2D)*(nu/nu_0)^beta
;	but that may make N_H unmanageably big, & conversion to those units hard
; CALLING:
;	gb = greybody( x, fdat, BETA=, KAPPA_0=, ... see kwarg list)
; INPUTS:
;	x = array of wavelengths or frequencies, in any units allowed for XUNIT.
;		defaults to wavelength in microns
;	fluxarr = 1D array of fluxes, or stack of images, at each wavelength
;		in either Jy / arcsec^2, MJy/sr, or W/m^2/hz
; KEYWORDS:
;	BETA = initial guess for emissivity
;	KAPPA_0 = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg*
;	NT = optional # of temp. components (only 1 working atm)
;	T_0 = scalar or length-2 array 
;	NCOL = initial guess for column density of dust, in either kg/m^2 or m^-2
;		default is kg/m^2
;	G2D = gas-to-dust ratio; default is 133
;	XUNIT = string units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = alphanumeric string units of fdat
;		accepts 'MJysr', 'Jyarcsec2', or 'Wm2hz' 
;		default is 'Jyarcsec2'
;	MASSCOL = Bool, 1 (True) if N_col is a mass density, 0 (False) if N_col is # density
;		default is 1 (#-density not working yet); if 0, G2D is set to its default
;	*note: 1 m^2/kg = 10 cm^2/g
; OUTPUTS:
;	Function returns array of Planck spectrum in units of Janskys/arcsec^2.
;	(Jansky = 1e-26 Watts/m^2/Hz)
; COMMON BLOCKS:
;	TBD
; EXTERNAL CALLS:
;	funit_conv
; HISTORY:
;	Written: Rebecca Pitts, 2016.
;---------------------------------------------------

function planck, nu, T
    	hc2 =  2 * double(!const.h) / double(!const.c)^2
    	hkt = double(!const.h) / ( double(!const.k) * T )
    	bb = hc2 * nu^3 / (exp( hkt * nu ) - 1)
	return, bb
  end

function mbb_opthin, nu, A
	;; A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = N_col, A[4] = G2D, A[5] = kappa_0
	;; assume mu = 2.8 & no need for omega
    	bb = planck(nu, A[0])
	plaw = ((nu / A[1])^A[2] )
	Sig = A[3] * double(2.8 * !const.mH) / A[4]
	return, [ bb * plaw * Sig * A[5]]
  end

;; use APEX SAM threshold to determine fg/bg contamination (maybe)
;;function rt1d_mbb, nudat, fldat,
;;	A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = omega, A[4] = N_col, A[5] = kappa_0
;;    	bb = planck(nu, A[0])
;;	return, [(nu / A[1])^A[2] ) * A[3] * A[4] * A[5] * bb]
;;  end

;;function PAHmodel,...

function x2nu, x, u
	case u of
	    'm': nudat = reverse((double(!const.c))/x)
	    'mm': nudat = reverse((double(!const.c) * 1e+3) / x)
	    'um': nudat = reverse((double(!const.c) * 1e+6) / x)
	    'nm': nudat = reverse((double(!const.c) * 1e+9) / x)
	    'Hz': nudat = double(x)
	    'MHz': nudat = double(x) * 1e+6
	    'GHz': nudat = double(x) * 1e+9
	    'THz': nudat = double(x) * 1e+12
	    else: message, 'Not a valid wavelength or frequency unit. Check spelling and case.'
          endcase
  	return, nudat
  end

function fucon, f, u, xu, invrs
    ;; doesn't matter what invrs is if not 0
    if (where(strmatch(xu,'*Hz',/fold_case)) eq -1) then rf = reverse(reform(double(f))) else $
     	  rf = reform(double(f)) ;;reverse(a,3) reverses array depthwise
    ;; if invrs=1, flux array should be already reversed so reversing again should set it right
    if (invrs eq 0) or ~keyword_set(invrs) then begin
	case u of
	    ;; 1 MJy/sr = 2.350443e-5 Jy/arcsec^2 = 10^-20 W/m^2/Hz
	    'Jyarcsec2': newf = (rf / double(2.350443e-5)) * 1e-20
	    'MJysr' : newf = rf * 1e-20
	    'Wm2Hz' : newf = rf
	    else: message, "Not a valid flux density unit. Enter 'MJysr', 'Jyarcsec2', or 'Wm2Hz'."
	  endcase
      endif else begin
	case u of
	    'Jyarcsec2': newf = (rf * double(2.350443e-5)) * 1e+20
	    'MJysr' : newf = rf * 1e+20
	    'Wm2Hz' : newf = rf
	    else: message, "Not a valid flux density unit. Enter 'MJysr', 'Jyarcsec2', or 'Wm2Hz'."
	  endcase
      endelse
  return, newf
  end

function greybody, x, fluxarr, XUNIT=xunit, FUNIT=funit, T_0=t_0, BETA=beta, $
			NCOL=ncol, MU=mu, KAPPA_0=kappa_0, G2D=g2d;;, $
			;;PAR_INIT=parmInit, PAR_MIN=Pmin, PAR_MAX=Pmax, PAR_FIT=parmFit, $
			;;NPARFIX=nparfix, PLOTFIT=plotfit, $
			;;INFO=info, VERBOSE=verbose, SHOW=show, MAXITER=maxit
			;;use NPARFIX to pick which parameters to fit

   ;;common IR_spectrum_fit, fitControl, WghtControl
   ;;common IR_spectrum_fit0, fitInitParams
   ;;common IR_spectrum_fit1, fitPlotOptions
   ;;common IR_spectrum_fit2, FitParLimits

	Nwav = N_elements(x)
	if (Nwav LE 0) then RETALL
	if (Nwav NE size(reform(fluxarr),/dimensions)) then begin
		message,"Error: single pixel stack extraction failed" $
			+ string(7b),/INFO
		RETALL,Nwav
          endif

    	if n_elements(xunit) lt 1 then xunit='um'
    	if n_elements(funit) lt 1 then funit='Jyarcsec2'
	nudat = x2nu(x,xunit)
    	fdat = fucon(fluxarr,funit,xunit,0)
	if funit eq 'Jyarcsec2' then omega = double(( !DTOR/3600 )^2) else omega = 1.0

	if n_elements(beta) ne 1 then begin
		if (n_elements(kappa_0) ne 1) then begin
			;;Default to Planck All-Sky values*
			;; *If D/G>100, kappa(250um) < 5.5 cm^2/g
			beta = 1.8
			kappa_0 = 0.55 ;; 0.55 m^2/kg from tau/Nh = 9.2e-26 cm^2
			nu0 = (double(!const.c) * 1e+6)/250
		  endif else beta = 2.0
	  endif
	if (beta eq 2.0) then begin
		if (n_elements(kappa_0) ne 1) then begin
			;;Second default to Herschel values
			;;(Magrini et al 2011, Davies et al 2012, Smith et al 2012)
			kappa_0 = 0.192
			nu0 = (double(!const.c) * 1e+6)/350
		  endif else begin
			if (kappa_0 gt 0.1) and (kappa_0 lt 0.3) then nu0 = (double(!const.c) * 1e+6)/350 else $
			if (kappa_0 gt 0.4) and (kappa_0 lt 0.6) then nu0 = (double(!const.c) * 1e+6)/250 else $
			STOP, 'opacity incompatible with built-in emissivity and fiducial frequency values.'
		  endelse	
	  endif
	;;lower beta should be better since that makes the power law curve wider & flatter

	;;Maybe one day I'll incorporate the ability to calculate the fiducial wavelength & kappa
	;; BUT IT IS NOT THIS DAY.

	if n_elements(ncol) ne 1 then ncol = 5.0e+21 * 1e+4 ;;H2mol/m^2, number density
	if n_elements(mu) ne 1 then mu = double(2.8) ;;mmw per unit H - assumes 71% H, 27% He, 2% Z
	;;from Arthur Cox 2000 reproduction of "Allen's Astrophysical Quantities" (1955) - everyone uses this
	;;alternative cited by Wikipedia doesn't seem to be reliable
	if n_elements(g2d) ne 1 then g2d = double(133.) ;;from Compiegne 2011
	if g2d lt 1 then g2d = 1/g2d
	
	if (n_elements(T_0) eq 1) OR (T_0[-1] eq 0) then begin ;if T_0 is INT, T_0[-1]=T_0[0]
	    Tc_0 = double(T_0[0])

	    ;; A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = omega, A[4] = N_col, A[5] = kappa_0
	    ;; 1st value in par_info to be replaced with Tc_0 later
	    par_info = replicate({value:0.D, fixed:1, limited:[1,1], limits:[0.D,0]}, 6)
	    par_info[*].value = [Tc_0,nu0,beta,ncol,g2d,kappa_0] ;; omega should always be fixed
	    par_info[0].fixed = 0
	    ;par_info[0].limits = [2.73,90]
	    ;par_info[1].limits = [min(nudat),max(nudat)]
	    ;par_info[2].limits = [0,5.0]
	    ;par_info[3].limits = [double(1.0e24),double(1.0e28)]
	    ;par_info[4].limits = [100.0,160.0]
	    ;par_info[5].limits = [0.01,1.0]

	    dstruct = {x:nudat,f:fdat}
	    pars = mpfit("mbb_opthin",functargs=dstruct,parinfo=par_info)
	    fmod = mbb_opthin(nudat, pars)
	  endif
	restruct = {nu:nudat, fdata:fdat, params:pars, fmodel:fmod}

return, restruct
end

;function greybody, x, imcube, XUNIT=xunit, FUNIT=funit, T_0=[20.,0], MU=mu$
;			BETA=beta, KAPPA_0=kappa_0, NCOL=ncol, G2D=g2d;;, $
;			;;PAR_INIT=parmInit, PAR_MIN=Pmin, PAR_MAX=Pmax, PAR_FIT=parmFit, $
;			;;NPARFIX=nparfix, PLOTFIT=plotfit, $
;			;;INFO=info, VERBOSE=verbose, SHOW=show, MAXITER=maxit
;			;;use NPARFIX to pick which parameters to fit
;
;   ;;common IR_spectrum_fit, fitControl, WghtControl
;   ;;common IR_spectrum_fit0, fitInitParams
;   ;;common IR_spectrum_fit1, fitPlotOptions
;   ;;common IR_spectrum_fit2, FitParLimits
;
;	;; probably need to move next 2 blocks to another function to loop over this one
;	Nwav = N_elements(x)
;	imd = size(imcube,/dimensions) ;;imd[-1] will be the same whether 1D or 3D - I checked
;	if (Nwav LE 0) then STOP
;
;	if (Nwav NE imd[-1]) then begin
;		message,"Error: Image cube depth =/= size of wavelength array" $
;			+ string(7b),/INFO
;		STOP,Nwav,imd
;          endif
;	
;return	
;end
	
