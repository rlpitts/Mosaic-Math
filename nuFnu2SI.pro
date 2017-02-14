function nuFnu2SI, x, f, xu, fu, BEAM=beam, XOFF=xoff
;+
; PURPOSE:
;	Convert any paired arrays of fluxes & frequencies/wavelengths to any
;	other units, using SI units as intermediaries
; ARGS:
;	X = array of wavelengths or frequencies. Note that
;	    program assumes you want them in increasing order
;	F = array of fluxes, either 1D or 3D, with the [last]
;	    dimension equal to the length of X
;	XU = alphanumeric string input unit of X (remove math symbols)
;	FU = string input unit of F
; KWARGS:
;	BEAMS = array of beam FWHMs for each frequency/wavelength
;		*required if either FU or FUOUT is 'Jy/beam'
;	XUOUT = output units of X; defaults to Hz
;	FUOUT = output units of F; defaults to W/m^2/Hz
; CALLING:
;	SIdata = nuFnu2SI(x, F, xu, fu)
;	*Assumes !const is defined & is in SI units
; OUTPUTS:
;	x & F in either SI units or units specified by iUOUT.
;	Recommend copying both inputs to new variables that
;	better reflect the output & don't override input
; EXTERNAL CALLS FROM:
;	greybody.pro
;	gbpstruct.pro
; HISTORY:
;	Written: Rebecca Pitts, 2017.
;---------------------------------------------------
  IF ~keyword_set(xoff) THEN BEGIN
    CASE xu OF
	'm': nudat = reverse((double(!const.c))/x)
	'mm': nudat = reverse((double(!const.c) * 1e+3) / x)
	'um': nudat = reverse((double(!const.c) * 1e+6) / x)
	'nm': nudat = reverse((double(!const.c) * 1e+9) / x)
	'Hz': nudat = double(x)
	'kHz': nudat = double(x) * 1e+3
	'MHz': nudat = double(x) * 1e+6
	'GHz': nudat = double(x) * 1e+9
	'THz': nudat = double(x) * 1e+12
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
 
   CASE 1 OF
	;; 1 MJy/sr = 2.350443e-5 Jy/arcsec^2 = 10^-20 W/m^2/Hz
	(where(strmatch(['Jyarcsec2','Jyasec2','Jy/asec2','Jy/arcsec2'],fu,/fold_case)) NE -1): newf = (rf / double(2.350443e-5)) * 1e-20
	(where(strmatch(['MJysr','MJy/sr'],fu,/fold_case)) NE -1) : newf = rf * 1e-20
	(where(strmatch(['Wm2Hz','W/m2/Hz'],fu,/fold_case)) NE -1) : newf = rf
        (where(strmatch(['Jy/beam','Jybeam'],fu,/fold_case)) NE -1) : BEGIN
           IF (n_elements(beam) EQ 2) THEN ba = double(!const.pi) * beam[0] * beam[1] / (4*alog(2)) ;;arcsec^2
           IF (n_elements(beam) EQ 1) THEN ba = double(!const.pi) * (beam^2) / (4*alog(2)) ;;arcsec^2
           IF ((n_elements(beam) LE 0) OR (n_elements(beam) GT 2)) THEN MESSAGE,$
		"ERROR: converting Jy/beam to W/m^2/Hz requires beam FWHM or [BMAJ,BMIN] in arcsec"
           newf = (rf / ba) * (1e-20 / double(2.350443e-5))
        END
	ELSE: message, "Not a valid flux density unit. Enter 'MJy/sr', 'Jy/asec2', 'Jy/beam', or 'W/m2/Hz'."
   ENDCASE

   IF ~keyword_set(xoff) THEN SIdata = {nu:nudat,flux:newf} ELSE SIdata = newf
    
   RETURN, SIdata
 END

    
