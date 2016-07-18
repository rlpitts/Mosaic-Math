;+
; NAME:
;	dust_Tmin
; PURPOSE:
;	Solve for the minimum Temperature of dust emitting over a distribution
;	of temperatures, thus the emission spectrum is more realistic.
;	The distribution of temperature is composed of two power laws
;	(for low and high temps. resp.) and depends on dust emissivity indices
;	(see function in file dust_spectrum.pro). The dust emission spectrum
;	of solution is stored in: common dust_Tmin_spec, dtd_spec.
; CALLING:
;	Tmin = dust_Tmin( Lumin, Tdust, Tcut, Tmax, $
;			EM=emiss_indices, W=wave, K=Kemit )
; INPUTS:
;	Lumin = the total Luminosity absorbed by the dust.
;	Tdust = temperature obtained by assuming all dust emission is at
;		single temperature (from function Dust_Temp).
;	Tcut = temperature at which dust emissivity index changes
;		(e.g. from 2 to 1 for graphite or 2 to 0 for silicates).
;	Tmax = maximum allowed dust temperature of distribution.
; KEYWORDS:
;	EMISS_INDEX = one or two element array of emissivity power law indices.
;	DENS_INDEX = optional, radial power law index of dust density variation,
;		e.g. DENS_INDEX=1 implies density goes as 1/r, default = 0.
;	LUMIN_INDEX = radial power law index of absorbed Luminosity, default=2.
;		Note: the negative of all above power law indices are used,
;			so normally specify indices as positive values.
;	KEMIT = dust emissivity (absorption cross-section).
;	WAVE = wavelenths in microns corresponding to KEMIT.
;	TOLERANCE = numerical tolerance for solution, default = 0.1.
; OUTPUTS:
;	Function returns minimum dust temperature of power-law distribution.
; COMMON BLOCKS:
;	common dust_Tmin_spec, dtd_spec		;resulting emission spectrum.
;	common dust_Tmin, Lum, Tx, Tm, dindex, wvdang
; EXTERNAL CALLS:
;	function Zbrent		(finds zero of a function)
;	function Trapez
;	function Dust_Spectrum
;	function dust_Tmin_spec	(included in this file)
; PROCEDURE:
;	Load common blocks, initialize function dust_spectrum, and solve
;	by matching absorbed and emitted Luminosities using function Zbrent.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function dust_Tmin_spec, Tmin

   common dust_Tmin, Lum, Tx, Tm, emindex, dindex, Lindex, wvdang
   common dust_Tmin_spec, dtd_spec

	dtd_spec = Dust_Spectrum( Tmin, Tx, Tm, E=emindex, D=dindex, L=Lindex )
	if !DEBUG then message, strconcat([Lum,Tmin,emindex,dindex,Lindex]),/INF

return, Lum - Trapez( dtd_spec, wvdang )
end
;-------------------------------------------------------------------------------
function dust_Tmin, Lumin, Tdust, Tcut, Tmax, TOLERANCE=Tol, DENS_INDEX=dinx, $
		EMISS_INDEX=emissinx, LUMIN_INDEX=Linx, WAVE=wave, KEMIT=Kemit

   common dust_Tmin, Lum, Tx, Tm, emindex, dindex, Lindex, wvdang

	if N_elements( Tol ) NE 1 then Tol = 0.1
	if N_elements( dinx ) NE 1 then dindex = 0 else dindex = dinx
	if N_elements( Linx ) NE 1 then Lindex = 2 else Lindex = Linx
	Lum = Lumin
	Tx = Tcut
	Tm = Tmax
	emindex = emissinx
	wvdang = wave * 1e4

	sp = Dust_Spectrum( Tdust, Tcut, Tmax, W=wave, K=Kemit, /INIT )

return, Zbrent( Tdust/(2.+(dindex>0)^4), Tdust, FUNC="dust_Tmin_spec", TOL=Tol )
end
