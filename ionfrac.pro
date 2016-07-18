;+
; NAME:
;	ionfrac
; PURPOSE:
;	Compute local ionization equilibrium, returning what fraction
;	of hydrogen is ionized, given density of Lyman continuum photons.
; CALLING:
;	fion = ionfrac( Dens_Lyc, DENS_H= )
; INPUTS:
;	Dens_Lyc = density of Lyman continuum photons (#/cm^3),
;		if a scalar, then we assume photons are all at Lyman edge.
; KEYWORDS:
;	DENS_H = density of hydrogen (default = 1/cm^3).
;	AION = ionization cross-section (default = 6.3e-18).
;	AREC = recombination coefficient (default = 3e-13).
;	FLUX_LYC = alternative way of specifying flux of Lyman continuum photons
;			(#/cm^2/sec) instead of the density of photons.
; OUTPUTS:
;	Function returns fraction of ionized gas.
; PROCEDURE:
;	Solve quadratic equation of local ionization balance.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function ionfrac, Dens_Lyc, DENS_H=Dens_H, SED=sed, WAVELENGTH=wavelens, $
				AION=aion, AREC=arec, NFLUX_LYC=nflux

	if N_elements( Dens_Lyc ) gt 0 then nflux = !cv.c * Dens_Lyc/(4*!PI)

	if N_elements( wavelens ) gt 1 then begin
		if N_elements( sed ) gt 1 then $
				nflux = sed * wavelens*1e-4/( !cv.h * !cv.c )
		if N_elements( nflux ) gt 1 then begin
			wang = wavelens*1e4
			piond = !cv.Lsun* Trapez( ismeuv( wang,1 )*nflux, wang )
		   endif
	 endif else begin
		if N_elements( aion ) ne 1 then aion = 6.3e-18
		piond = aion * double( nflux )
	  endelse

	if N_elements( piond ) LE 0 then begin
		print,"syntax:  fion = ionfrac( Dens_Lyc, DENS_H= )"
		print,"    or:  fion = ionfrac( NFLUX_LYC=, DENS_H= )"
		print,"    or:  fion = ionfrac( SED=, WAVE=, DENS_H= )"
		return,0
	   endif

	if N_elements( Dens_H ) LE 0 then begin
		Dens_H = 1
		message,"assuming n(H) = 1",/INFO
	   endif
help,piond

	piond = piond/Dens_H
	if N_elements( arec ) ne 1 then arec = 3e-13

return, ( sqrt( 4*arec/piond + 1 ) - 1 ) * piond/(2*arec)
end
