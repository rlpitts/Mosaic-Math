;+
; NAME:
;	Rad_Stromgren
; PURPOSE:
;	Compute Stromgren radius of an idealized spherical HII region
;	with central source of ionizing photons.
;
; CALLING:
;	radius = Rad_Stromgren( Log_N_Lyc, DENSITY_H = )
;
; INPUTS:
;	Log_N_Lyc = Log base 10 of number of Lyman continuum ionizing photons.
;
; KEYWORDS:
;	DENSITY_H = density of Hydrogen atoms per cubic cm,  default=1.
;
;	RECOMBINATION = recombination rate coefficient, default = 3e-13.
;
; OUTPUTS:
;	Function returns Stromgren radius of HII region in parsecs.
;
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function Rad_Stromgren, Log_N_Lyc, RECOMBINATION=arec, DENSITY_H=Dens_H, $
					SED=sed, WAVELENGTH=wavelens

	if N_elements( sed ) gt 1 then begin
		wlyc = where( wavelens LT 0.0912, nw )
		if( nw LE 1 ) then begin
			Log_N_Lyc = 0
			return,0
		   endif
		hc = !cv.h * !cv.c * 1e4	;units = erg-microns
		hfinv = wavelens(wlyc)/hc	; = 1/(h*freq) = ergs^(-1)
		Log_N_Lyc = alog10( Trapez( sed * hfinv, wavelens ) ) + 4 $
				+ alog10( !cv.Lsun )
	   endif

	if N_elements( arec ) NE 1 then arec = 3e-13
	if N_elements( Dens_H ) NE 1 then Dens_H = 1

return,10^(( Log_N_Lyc + alog10(3) - alog10(4*!PI*arec) - 2*alog10(Dens_H) )/3 $
		- alog10(3.086) - 18 )
end
