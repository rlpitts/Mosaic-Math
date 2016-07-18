;+
; NAME:
;	dTemp_Prob
; PURPOSE:
;	Compute the probabilities of given dust temperatures around a star,
;	assuming a distribution of temperatures of power law form.
;	Normalization is based on the range (min-max) of input temperatures.
;	Based on optically thin theory of Luminosity vs. radius and
;	temperature vs. absorbed Luminosity, coupled with density vs. radius.
; CALLING:
;	probability = dTemp_Prob( Tgrid )
; INPUTS:
;	Tgrid = array of dust temperatures of interest, degrees Kelvin.
;		Note that Tmin and Tmax are obtained from this grid.
; KEYWORDS:
;	EMISS_INDEX = one or two element array of emissivity power law indices.
;		e.g. [2,1] for graphite or [2,0] for silicates, default = [2,1].
;	TCUT = temperature at which dust emissivity index changes
;		if 2 indices are specified, default = 100 K.
;	DENS_INDEX = radial power law index of dust density variation,
;		e.g. DENS_INDEX=1 implies density goes as 1/r, default = 0.
;	LUMIN_INDEX = radial power law index of absorbed Luminosity, default=2.
;		Note: the negative of all above power law indices are used,
;			so normally specify indices as positive values.
;	/AVERAGE : return instead the scalar average temperature of distrib.
;		In this case you should pass Tgrid = [ Tmin, Tmax ].
; OUTPUTS:
;	Function returns probabilities of input dust temperatures,
;	normalized over the given range of temperatures.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function dTemp_Prob, Tgrid, TCUT=Tcut, EMISS_INDEX=emindex, $
			LUMIN_INDEX=Lumindex, DENS_INDEX=denindex, AVERAGE_T=avt

	if N_elements( Tcut ) NE 1 then Tcut = 100
	if N_elements( emindex ) LE 0 then emindex = [ 2, 1 ]
	if N_elements( Lumindex ) LE 0 then Lumindex = 2
	if N_elements( denindex ) LE 0 then denindex = 0

	gama = 1 + float( 3 - denindex )*( 4 + emindex )/Lumindex
	if !DEBUG then $
		message,"T power Law(s) = " + string( -gama, F="(2F8.2)" ),/INFO

	if N_elements( Tgrid ) LE 0 then begin
		print,"syntax:	p_T = dTemp_Prob( Tgrid )"
		print,"keywords:		TCUT="
		print,"			EMISS_INDEX="
		print,"			LUMIN_INDEX="
		print,"			DENS_INDEX="
		return,0
	  endif else Tmin = min( Tgrid, MAX=Tmax )

	if N_elements( gama ) EQ 2 then begin

		if (Tmin GE Tcut) then $
			return, dTemp_Prob( Tgrid, EM=emindex(1), $
					LUM=Lumindex, DEN=denindex, AV=avt )
		if (Tmax LE Tcut) then $
			return, dTemp_Prob( Tgrid, EM=emindex(0), $
					LUM=Lumindex, DEN=denindex, AV=avt )
		p = 1 - gama
		T1 = [ Tcut, Tmax ]
		T2 = [ Tmin, Tcut ]
		tc = ( double( T1 )^p - double( T2 )^p )/p
		tf = Tcut^(gama(1)-gama(0))
		tc(1) = tc(1) * tf
		ac = 1/total( tc )
		ac = [ ac, ac * tf ]
		if !DEBUG then print,tc,tf,ac
		if !DEBUG GT 1 then stop

		if keyword_set( avt ) then begin
			p = 2 - gama
			bc = ac/p
			return, total( bc * ( double( T1 )^p - double( T2 )^p ))
		 endif else begin
			pt = dblarr( N_elements( Tgrid ) )
			w = where( Tgrid LE Tcut, nw )
			if (nw GT 0) then pt(w) = ac(0) * Tgrid(w)^(-gama(0))
			w = where( Tgrid GT Tcut, nw )
			if (nw GT 0) then pt(w) = ac(1) * Tgrid(w)^(-gama(1))
			return,pt
		  endelse

	 endif else begin

		gama = gama(0)
		p = 1 - gama
		ac = p/( double( Tmax )^p - double( Tmin )^p )
		if !DEBUG then print,ac
		if !DEBUG GT 1 then stop

		if keyword_set( avt ) then begin
			p = 2 - gama
			return, (ac/p) * ( double( Tmax )^p - double( Tmin )^p )
		 endif else return, ac * Tgrid^(-gama)
	  endelse
end
