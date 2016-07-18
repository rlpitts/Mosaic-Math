;+
; NAME:
;	dust_Temp
; PURPOSE:
;	Compute the radiative equilibrium temperature of dust by matching
;	absorbed luminosity and emission integrated over given wavelengths.
;	The emision spectrum is placed in common dust_emit_spec for retrieval.
; CALLING:
;	Temperature = dust_Temp( Lum_abs, wavel, kabs )
; INPUTS:
;	Lum_abs = luminosity absorbed by the dust, in ergs.
;	wavel = wavelengths in microns.
;	kabs = absorption cross-section at wavelengths (=emissivity).
; KEYWORDS:
;	RADIUS_SPHERE = radius (in cm) of homogenous sphere
;		in which given Luminosity is absorbed,
;		and then units are assumed to be solar luminosities,
;		and kabs is assumed to be the radial optical depth of sphere.
;		(this keyword is optional, just for a specific application)
;	TOLERANCE = numerical tolerance for solution, default = 0.1.
;	MINT = mininum T, default = 3.
;	MAXT = maximum T, default = 2000.
; OUTPUTS:
;	Function returns effective dust temperature in degrees kelvin.
; COMMON BLOCKS:
;	common dust_emission, sigma, wavelang	;for internal purposes.
;	common dust_emit_zero, Lum_emit
;	common dust_emit_spec, dust_spec	;the resultant spectrum
; EXTERNAL CALLS:
;	function Zbrent		(finds zero of a function)
;	function Planck
;	function Trapow
;	function dust_emit_zero	(included in this file)
; PROCEDURE:
;	Find zero of:  Lum_abs - 4 * Integral( kabs * Planck( Temperature ) )
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function dust_emit_zero, Temperature

   common dust_emission, sigma, wavelang
   common dust_emit_zero, Lum_emit
   common dust_emit_spec, dust_spec

	dust_spec = 4 * sigma * Planck( wavelang, Temperature )
	if !DEBUG then message,strconcat( [Lum_emit, Temperature] ),/INFO

return, Lum_emit - float( Trapow( dust_spec > 1e-37, wavelang ) )
end

;-----------------------------------------------------------------------------

function dust_Temp, Lum_abs, wavel, kabs, RADIUS_SPHERE=rs, TOLERANCE=Tol, $
						MINT=minT, MAXT=maxT
   common dust_emission, sigma, wavelang
   common dust_emit_zero, Lum_emit

	if keyword_set( rs ) then $
		sigma = kabs * (rs/3.826e33) * (4*!pi*rs/3) $
	  else	sigma = kabs

	Lum_emit = Lum_abs
	wavelang = 1e4 * wavel
	if N_elements( Tol ) NE 1 then Tol = 0.1
	if N_elements( minT ) NE 1 then minT = 3.0
	if N_elements( maxT ) NE 1 then maxT = 2000.0

return, Zbrent( minT, maxT, FUNC="dust_emit_zero", TOL=Tol )
end
