;+
; NAME:
;	greybody.pro
; PURPOSE:
;	Compute greybody emission spectrum (Planck function)
;	in units of Janskys/arcsec^2 at for an array of wavelengths.
;	Uses the following relationships:
;		kappa = kappa_0*(wl_0/wl)^beta
;		kappa_0*N_H = tau_0, where N = col. density & tau << 1
;	If tau = absorption cross-section & N_H = # of H atoms, then
;		kappa*m_H*(dust/gas) = tau/N_H = (tau_0/N_H)*(wl_0/wl)^beta
;	but that may make N_H unmanageably big, & conversion to those units hard
; CALLING:
;	gb = greybody( x, fdat, BETA=, KAPPA_0=, NCOL=)
; INPUTS:
;	x = array of wavelengths, in microns.
;
;	fdat = 1D array of fluxes, or stack of images, at each wavelength, in Jy / arcsec^2.
; KEYWORDS:
;	BDIMS = scalar or length-2 array describing beam size, in arcsec
;	BETA = initial guess for emissivity
;	KAPPA_0 = initial guess for opacity at fiducial frequency nu_0 (tbd), in m^2/kg
;	NCOL = initial guess for column density of dust, in kg/m^2
;	XUNIT = string specifying units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = alphanumeric string, either 'MJysr' or 'Jyarcsec2', specifying units of fdat
;		default is 'Jyarcsec2'
;	T0 = optional initial guess for blackbody peak temperature; default is 20 K
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

function mplanck, nu, A
	;; A[0] = T, A[1] = nu_0, A[2] = beta, A[3] = omega, A[4] = N_col, A[5] = kappa_0, 
	;; 	A[6] = coeffs of Planck function
	xT = exp( double(!const.h) * nu / ( double(!const.k) * A[0] ) ) - 1
	return, ( (nu / A[1])^A[2] ) * A[3] * A[4] * A[5] * A[6] / xT

function rt1d_mbb, nudat, fldat,
    hc2 =  2 * double(!const.h) / double(!const.c)^2
    hkt = double(!const.h) / ( double(!const.k) * T0 )
    bb = hc2 * nudat^3 / (exp( hkt * nudat ) - 1)
    return, '...'

function PAHmodel


function greybody, x, fdat, BDIMS=bdims, XUNIT=xunit, FUNIT=funit, TC_0=Tc_0, $
			BETA=beta, KAPPA_0=kappa_0, NCOL=ncol, D2G=d2g;;, $
			;;TH_0=Th_0, NPARFIX=nparfix, PLOTFIT=plotfit, $
			;;INFO=info, VERBOSE=verbose, SHOW=show, MAXITER=maxit
			;;use NPARFIX to pick which parameters to fit

   ;;common IR_spectrum_fit, fitControl, WghtControl
   ;;common IR_spectrum_fit0, fitInitParams
   ;;common IR_spectrum_fit1, fitPlotOptions
   ;;common IR_spectrum_fit2, FitParLimits

	Nwav = N_elements(x)
	imd = size(fdat,/dimensions) ;;imd[-1] will be the same whether 1D or 3D - I checked
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
    
    	if (xunit eq 'm') || (xunit eq 'mm') || (xunit eq 'um') || (xunit eq 'nm') then rf = reverse(reform(double(fdat)),3) else $
     	 rf = reform(double(fdat))
	;;reverse(a,3) reverses array depthwise

    	case funit of
	    ;; 1 MJy/sr = 2.350443e-5 Jy/arcsec^2 = 10^-20 W/m^2/Hz
	    'Jyarcsec2': fdat = (rf / double(2.350443e-5)) * 1e-20
	    'MJysr' : fdat = rf * 1e-20
	    else: message, "Not a valid flux density unit. Enter 'MJysr' or 'Jyarcsec2'."
	 endcase
	;;want everything in order of increasing frequency

	sqas2rad = ( !DTOR/3600 )^2
	if n_elements(bdims) lt 1 then omega = 1.0 * sqas2rad else $ ;;use this when units are arcsec/px
	if n_elements(bdims) eq 1 then omega = double(!const.pi) * sqas2rad * bdims^2/(4 * alog(2)) else $
	if n_elements(bdims) eq 2 then omega = double(!const.pi) * sqas2rad * bdims[0] * bdims[1]/(4 * alog(2)) else $
	 print,'Error: enter beam dimensions in arcsec.' ;;sq arcsec converstion already

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
			kappa_0 = 1.92
			nu0 = (double(!const.c) * 1e+6)/350
		  endif else begin
			if (kappa_0 gt 0.1) and (kappa_0 lt 0.3) then nu0 = (double(!const.c) * 1e+6)/350 else $
			if (kappa_0 gt 0.4) and (kappa_0 lt 0.6) then nu0 = (double(!const.c) * 1e+6)/250 else $
			 print, 'opacity incompatible with built-in emissivity and fiducial frequency values.'
			return,0
		  endelse	
	  endif
	;;lower beta should be better since that makes the power law curve wider & flatter

	;;Maybe one day I'll incorporate the ability to calculate the fiducial wavelength
	;; BUT IT IS NOT THIS DAY.
	
	if n_elements(ncol) ne 1 then ncol = (5e+21) * 1e+4 ;;H2mol/m^2

	;;coeff = (2 * !const.h / !const.c^2)
	;;mpfit or lmfit?

return

end
	
