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
;	gb = greybody( x, fdat, BETA=, KAPPA_0=, T0=, NCOL=, G2D=, XUNIT=, FUNIT=, NUNIT=)
; INPUTS:
;	x = array of wavelengths or frequencies, in any units allowed for XUNIT.
;		defaults to wavelength in microns
;	fdat = 1D array of fluxes, or stack of images, at each wavelength, in Jy / arcsec^2.
; KEYWORDS:
;	BETA = initial guess for emissivity
;	KAPPA_0 = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg*
;	NT = optional # of temp. components (only 1 working atm)
;	NCOL = initial guess for column density of dust, in either kg/m^2 or m^-2
;		default is kg/m^2
;	G2D = gas-to-dust ratio; default is 133
;	XUNIT = string specifying units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = alphanumeric string, either 'MJysr' or 'Jyarcsec2', specifying units of fdat
;		default is 'Jyarcsec2'
;	NUNIT = string, either 'mass' or 'num', indicating if N_col is a mass or # density
;		if 'num' and G2D isn't specified, G2D is set to its default
;		default is mass-density; #-density not working yet
;	*note: 1 m^2/kg = 10 cm^2/g
; OUTPUTS:
;	Function returns array of Planck spectrum in units of Janskys/arcsec^2.
;	(Jansky = 1e-26 Watts/m^2/Hz)
; COMMON BLOCKS:
;	TBD
; EXTERNAL CALLS:
;	TBD (may use this to make function to call it in a loop)
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
	;; A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = omega, A[4] = N_col, A[5] = kappa_0
    	bb = planck(nu, A[0])
	return, [((nu / A[1])^A[2] ) * A[3] * A[4] * A[5] * bb]
  end

;;function rt1d_mbb, nudat, fldat,
;;	A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = omega, A[4] = N_col, A[5] = kappa_0
;;    	bb = planck(nu, A[0])
;;	return, [(nu / A[1])^A[2] ) * A[3] * A[4] * A[5] * bb]
;;  end

;;function PAHmodel,...


function greybody, x, fluxarr, XUNIT=xunit, FUNIT=funit, NUNIT=nunit, TCOMPS=tcomps, $
			BETA=beta, KAPPA_0=kappa_0, NCOL=ncol, G2D=g2d;;, $
			;;PAR_INIT=parmInit, PAR_MIN=Pmin, PAR_MAX=Pmax, PAR_FIT=parmFit, $
			;;TC_0=Tc_0, TH_0=Th_0, NPARFIX=nparfix, PLOTFIT=plotfit, $
			;;INFO=info, VERBOSE=verbose, SHOW=show, MAXITER=maxit
			;;use NPARFIX to pick which parameters to fit

   ;;common IR_spectrum_fit, fitControl, WghtControl
   ;;common IR_spectrum_fit0, fitInitParams
   ;;common IR_spectrum_fit1, fitPlotOptions
   ;;common IR_spectrum_fit2, FitParLimits

	Nwav = N_elements(x)
	imd = size(fluxarr,/dimensions) ;;imd[-1] will be the same whether 1D or 3D - I checked
	if (Nwav LE 0) then return,0

	if (Nwav NE imd[-1]) then begin
		message,"Error: Image cube depth =/= size of wavelength array" $
			+ string(7b),/INFO
		print,Nwav,imd
		return,0
          endif

    	if n_elements(xunit) lt 1 then xunit='um'
    	if n_elements(funit) lt 1 then funit='Jyarcsec2'

    	case xunit of
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
    
    	if (xunit eq 'm') || (xunit eq 'mm') || (xunit eq 'um') || (xunit eq 'nm') then rf = reverse(reform(double(fluxarr)),3) else $
     	  rf = reform(double(fluxarr))
	;;reverse(a,3) reverses array depthwise

    	case funit of
	    ;; 1 MJy/sr = 2.350443e-5 Jy/arcsec^2 = 10^-20 W/m^2/Hz
	    'Jyarcsec2': fdat = (rf / double(2.350443e-5)) * 1e-20
	    'MJysr' : fdat = rf * 1e-20
	    'Wm2Hz' : fdat = rf
	    else: message, "Not a valid flux density unit. Enter 'MJysr', 'Jyarcsec2', or 'Wm2Hz'."
	  endcase

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
			 print, 'opacity incompatible with built-in emissivity and fiducial frequency values.'
			return,0
		  endelse	
	  endif
	;;lower beta should be better since that makes the power law curve wider & flatter

	;;Maybe one day I'll incorporate the ability to calculate the fiducial wavelength & kappa
	;; BUT IT IS NOT THIS DAY.

	if n_elements(ncol) ne 1 then ncol = 5.0e+21 * 1e+4 ;;H2mol/m^2
	if n_elements(nunit) ne 1 then nunit = 'mass'
	if nunit eq 'mass' then begin
	    n_md = double(ncol)
	  endif else begin
	    if nunit eq 'num' then begin
		mu = 2.8 ;;mmw per unit H - assumes 71% H, 27% He, 2% everything else
		mh = double(1.673534e-27)
		if n_elements(g2d) ne 1 then g2d = 133
		n_md = mu*mh*double(ncol)/g2d
	      endif
	  endelse
	
	;; A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = omega, A[4] = N_col, A[5] = kappa_0
	;; 1st value in par_info to be replaced with Tc_0 later
	par_info = replicate({value:0.D, fixed:1, limited:[1,1], limits:[0.D,0]}, 6)
	par_info[*].value = [20,nu0,beta,omega,ncol,kappa_0] ;; omega should always be fixed
	par_info[0].fixed = 0
	par_info[0].limits = [2.73,90]
	par_info[1].limits = [min(nudat),max(nudat)]
	par_info[2].limits = [0,5.0]
	par_info[4].limits = [double(1.0e24),double(1.0e28)]
	par_info[5].limits = [0.01,1.0]

	if (n_elements(tcomps) ne 1) or (tcomps eq 1) then begin
	    dstruct = {x:nudat,f:fdat}
	    pars = mpfit("mbb_opthin",functargs=dstruct,parinfo=par_info)
	    fmodel = mbb_opthin(nudat, pars)
	  endif

return, [nudat, fmodel]
end
	
