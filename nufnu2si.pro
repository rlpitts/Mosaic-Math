function nufnu2si, x, f, xu, fu, PXSZ=pxsz, BEAM=beam, XOFF=xoff
;+
; PURPOSE:
;	Convert any paired arrays of fluxes & frequencies/wavelengths to 
;	Hz and W/(m^2*Hz*sr)
; ARGS:
;	X = array of wavelengths or frequencies. Note that
;	    program assumes you want them in increasing order
;	F = array of fluxes, either 1D or 3D, with the [last]
;	    dimension equal to the length of X
;	XU = alphanumeric string input unit of X (remove math symbols)
;	FU = string input unit of F
; KWARGS:
;	BEAM = array of beam FWHMs for each frequency/wavelength
;		*required if FU is 'Jy/beam'
;	XOFF = Boole; if set, return flux densities alone
; CALLING:
;	SIdata = nufnu2si(x, F, xu, fu)
;	*Assumes !const is defined & is in SI units
; OUTPUTS:
;	Recommend copying both inputs to new variables that
;	better reflect the output & don't override input
; EXTERNAL CALLS FROM:
;	greybody.pro
;	iter_gb3d.pro
; HISTORY:
;	Written: Rebecca Pitts, 2017.
;---------------------------------------------------
  IF ~keyword_set(xoff) THEN BEGIN
    CASE xu OF
	'm': nudat = reverse((double(!const.c))/x)
	'mm': nudat = reverse((double(!const.c) * 1e3) / x)
	'um': nudat = reverse((double(!const.c) * 1e6) / x)
	'nm': nudat = reverse((double(!const.c) * 1e9) / x)
	'Hz': nudat = double(x)
	'kHz': nudat = double(x) * 1e3
	'MHz': nudat = double(x) * 1e6
	'GHz': nudat = double(x) * 1e9
	'THz': nudat = double(x) * 1e12
	ELSE: MESSAGE, 'Not a valid wavelength or frequency unit. Check spelling and case.'
     ENDCASE
   ENDIF

   ndims = size(f,/n_dimensions)
   shape = size(f)
   IF ((n_elements(f) GT 1) && (where(strmatch(xu,'*Hz',/fold_case)) EQ -1)) THEN BEGIN
      IF ((ndims GT 1) AND (shape[1] GT 1) AND (shape[1] GT 1)) THEN rf = reverse(reform(double(f)),ndims) $
      ELSE rf = reverse(reform(double(f)))
   ENDIF ELSE BEGIN
      rf = reform(double(f))
   ENDELSE

   IF (N_elements(pxsz) NE 0) THEN BEGIN
	IF (N_elements(pxsz) EQ 1) THEN pxa = pxsz^2 ELSE pxa = pxsz[0]*pxsz[1]
	;;[asec^2/px] assumes px dimensions are already in asec
   ENDIF
   sr2asec = ( double(!const.pi) / (180.d*3600.d) )^2
 
   CASE 1 OF
	;; 1 MJy/sr = 2.350443e-5 Jy/arcsec^2 = 10^-20 W/m^2/Hz/sr
        ;; 1 Jy = 10^-26 W/m^2/Hz
	(where(strmatch(['Jyarcsec2','Jyasec2','Jy/asec2','Jy/arcsec2'],fu,/fold_case)) NE -1): newf = rf * (1e-26 / sr2asec)
	(where(strmatch(['MJysr','MJy/sr'],fu,/fold_case)) NE -1): newf = rf * 1e-20
	(where(strmatch(['Wm2Hz','W/m2/Hz','Wm2Hzsr','W/m2/Hz/sr','W/m2Hzsr'],fu,/fold_case)) NE -1): newf = rf
        (where(strmatch(['Jy/beam','Jybeam'],fu,/fold_case)) NE -1): BEGIN
           IF (n_elements(beam) EQ 2) THEN ba = double(!const.pi) * double(beam[0]) * double(beam[1]) / (4.*alog(2.)) ;;arcsec^2
           IF (n_elements(beam) EQ 1) THEN ba = double(!const.pi) * double(beam^2) / (4.*alog(2.)) ;;arcsec^2
           IF ((n_elements(beam) LE 0) OR (n_elements(beam) GT 2)) THEN MESSAGE,$
		"ERROR: converting Jy/beam to W/m^2/Hz/sr requires beam FWHM or [BMAJ,BMIN] in arcsec"
           newf = (rf / ba) * (1e-26 / sr2asec)
        END
	(where(strmatch(['Jypx','Jypixel','Jy/pixel','Jy/px'],fu,/fold_case)) NE -1): BEGIN
	   IF (N_elements(pxsz) EQ 0) THEN message, "Error: conversion from Jy/pixel requires pixel dimensions arcsec"
	   IF (N_elements(pxsz) EQ 1) THEN pxa = pxsz^2 ;;[asec^2/px] assumes px dimensions are already in asec
	   IF (N_elements(pxsz) EQ 2) THEN pxa = pxsz[0]*pxsz[1]
	   IF (N_elements(pxsz) GT 2) THEN message, "Error: voxels and higher dimensions not handled"
	   newf = (rf / pxa) * (1e-26 / sr2asec)
	END
	ELSE: message, "Not a valid flux density unit. Enter 'MJy/sr', 'Jy/asec2', 'Jy/beam', or 'W/m2/Hz/sr'."
   ENDCASE

   IF ~keyword_set(xoff) THEN SIdata = {nu:nudat,flux:newf} ELSE SIdata = newf
    
   RETURN, SIdata
 END

    
