;+
; NAME:
;	Dust_Spectrum
; PURPOSE:
;	Compute spectrum of emission from dust heated by point source(s)
;	assuming a distribution of temperatures of power law form.
;	Based on optically thin theory of Luminosity vs. radius and
;	temperature vs. absorbed Luminosity, coupled with density vs. radius.
; CALLING:
;	spectrum = Dust_Spectrum( Tmin, Tcut, Tmax )
; INPUTS:
;	Tmin = minimun temperature of dust (default = 2.7 Kelvin).
;	Tcut = turnover temperature between 2 emissivity indices (default=100).
;	Tmax = maximum temperature of dust (default = 1500 Kelvin).
; KEYWORDS:
;	EMISS_INDEX = emissivity power law index or 2 indices, default = [2,1].
;	DENS_INDEX = optional, radial power law index of dust density variation,
;		e.g. DENS_INDEX=1 implies density goes as 1/r, default = 0.
;	LUMIN_INDEX = radial power law index of absorbed Luminosity, default=2.
;		Note: the negative of all above power law indices are used,
;			so normally specify indices as positive values.
;	KABS = absorption cross-section of dust vs. wavelength (=emissivity).
;	WAVELENS = wavelengths (microns) of dust cross-section.
;	BBMATRIX = optional matrix of Planck spectra for distribution of
;		temperatures (default is to compute and store in common).
;	TGRID = optional temperature grid for distribution integration,
;		or if scalar, # of elements in default temperature grid.
;	/INIT : just load the keyword parameters into common block,
;		setup/compute BB matrix, and return,0.
; OUTPUTS:
;	Function returns spectrum (ergs/Angstrom/sec) of dust emission.
; COMMON BLOCKS:
;	common Dust_Spectrum, Kemit, wvang, Tg, BBmat
; EXTERNAL CALLS:
;	function dTemp_Prob
;	function Planck
;	function Trapow
; PROCEDURE:
;	Temperature distribution is obtained from function dTemp_Prob.
;	Carefully integrate over temperature distribution of Planck spectra.
;	The usual factor of 4*!pi steradians for emission is just 4 since
;	!pi is already included in Kemit = Kabs.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function Dust_Spectrum, Tmin, Tcut, Tmax, EMISS_INDEX=emindex, INIT=init, $
			LUMIN_INDEX=Lumindex, DENS_INDEX=denindex, KABS=Kabs,$
			WAVELENS=wavelens, TGRID=Tgrid, BBMATRIX=BBmatrix

   common Dust_Spectrum, Kemit, wvang, Tg, BBmat

	if N_elements( wavelens ) GT 1 then wvang = wavelens * 1e4
	if N_elements( Kabs ) GT 1 then Kemit = Kabs

	if N_elements( Tgrid ) GT 1 then begin
		Tg = Tgrid
		if N_elements( BBmatrix ) GT 1 then BBmat= transpose( BBmatrix )
	  endif else begin
		if N_elements( Tgrid ) EQ 1 then nt = Tgrid(0) else nt = 300
		Tg = 10^( aLog10(2.7) + aLog10(2000/2.7) * findgen(nt)/(nt-1))
	   endelse

	ntemp = N_elements( Tg )
	nwave = N_elements( wvang )
	sz = size( Kemit )
	if (nwave NE sz(1)) then message,"check # of wavelengths",/INFO
	if sz(0) GT 1 then ncomp=sz(2) else ncomp=1

	sz = size( BBmat )
	if (sz(0) NE 2) OR (ntemp NE sz(1)) OR (nwave NE sz(2)) then begin
		message,"computing Planck function (BB) matrix",/INFO
		BBmat = fltarr( nwave, ntemp )
		for i=0,ntemp-1 do BBmat(*,i) = Planck( wvang, Tg(i) )
		BBmat= transpose( BBmat )
	   endif

	if keyword_set( init ) then return,0

	if N_elements( Tmin ) NE 1 then Tmin = 2.7
	Tmin = Tmin > 2.7
	m = min( abs( Tg - Tmin ), imin )	;substitute for min Temp.

	if (Tg(imin) NE Tmin) then begin
		Tg(imin) = Tmin
		BBmat(imin,*) = reform( Planck( wvang, Tmin ), 1, nwave )
	   endif

	if N_elements( Tmax ) NE 1 then Tmax = 1500
	m = min( abs( Tg - Tmax ), imax )	;substitute for max Temp.

	if (Tg(imax) NE Tmax) then begin
		Tg(imax) = Tmax
		BBmat(imax,*) = reform( Planck( wvang, Tmax ), 1, nwave )
	   endif

	if N_elements( Tcut ) NE 1 then Tcut = 100
	if (Tcut GT Tmin) AND (Tcut LT Tmax) then begin
		m = min( abs( Tg - Tcut ), icut )	;subs for cut Temp.
		if (Tg(icut) NE Tcut) then begin
			Tg(icut) = Tcut
			BBmat(icut,*) = reform( Planck( wvang, Tcut ), 1,nwave )
		   endif
	   endif

	spectrum = fltarr( nwave, ncomp )
	if !DEBUG then help,imin,imax
	Ts = Tg(imin:imax)
	pT = dTemp_Prob( Ts, EM=emindex, LUM=Lumindex, DEN=denindex, TCUT=Tcut )

	for i=0,ncomp-1 do begin
	  for k=0,nwave-1 do begin	;integrate over Temp Distribution:
	     bT = BBmat(imin:imax,k)
	     w = where( bT GT 0, nw )
	     if (nw GT 1) then begin
		s = w(0)
		e = w(nw-1)		;use a connected set for integral.
		bT = bT(s:e) > 1e-37
		spectrum(k,i) = 4 * Kemit(k,i) * Trapow( bT * pT(s:e), Ts(s:e) )
	      endif
	    endfor
	  endfor

return, spectrum
end
