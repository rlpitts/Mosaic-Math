;+
; NAME:
;	fractal_synth
; PURPOSE:
;	Synthesize a fractal embedded in either 2-D (a fractal curve)
;	or in 3-D (a fractal surface) using power-law spectral technique.
; CALLING:
;	function fractal_synth, Hpar, ft
; INPUTS:
;	Hpar = H parameter determining the theoretical fractal dimension D,
;		should normally be in the range 0 to 1 so that D = E - H,
;		where E = embedding dimension (2 or 3). The theoretical fractal
;		dimension is approached as number of frequencies goes infinite.
;		The power-law exponent of power spectral density = -(2*H+E-1).
;		Statistically, the average variance of points on the fractal
;		scales as the distance between the points to the power 2*H.
; KEYWORDS:
;	EDIM = embedding dimension E, default=3, giving fractal surface (image).
;	NFREQ = number of frequencies to use in spectral synthesis, default=64.
;	STD = optional st.dev. of random noise added to amplitudes, default=0.
;	PHASE = optional input/output, random phases between 0 and 2*!pi.
;		Keyword can be used to keep the same phases while changing H.
; OUTPUTS:
;	ft = optional, Fourier spectrum used to generate the synthesized image.
;
;		Function returns a vector of size = 2 * NFREQ,
;			or image of size = ( 2 * NFREQ, 2 * NFREQ ).
; EXAMPLES:
;	display 128x128 image of dimension 2.5 :
;		tvscl, fractal_synth()
;
;	display 512x512 image of dim = 2.3 :
;		tvscl, fractal_synth( 0.7, NF=256 )
;
;	plot 512 pnt. curve of dim = 1.9 (and save phases into p) :
;		plot, fractal_synth( 0.1, E=2, NF=256, PH=p )
;
;	overplot curves with range of fractal dimensions from 1 to 2 :
;		for i=0,10 do oplot, fractal_synth( i/10., E=2, NF=256, PH=p )
; 
; EXTERNAL CALLS:
;	function rfftinv	(inverse FFT giving a real-valued function).
; COMMON BLOCKS:
;	common fractal_synth, su,sn   (internal use, seeds of randomu & randomn)
; PROCEDURE:
;	Generate 1D/2D spectral amplitudes with power-law exponent -(2*H+E-1)/2
;	and distribute the phases uniformly from 0 to 2*!pi, then inverse FFT.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function fractal_synth, Hpar, ft, EDIM=edim, NFREQ=nfreq, STD=std, PHASE=phase

   common fractal_synth, su,sn

	if N_elements( edim ) NE 1 then edim = 3
	ndim = edim-1
	if N_elements( nfreq ) LE 0 then nfreq = 64
	if N_elements( nfreq ) NE ndim then nfreq = replicate( nfreq(0), ndim )
	if N_elements( Hpar ) NE 1 then Hpar = 0.5

	CASE edim OF

	  2: BEGIN
		message,"fractal dimension = " + $
			strtrim( ((edim-hpar)>1)<2, 2 ),/INFO
		power = -( Hpar + 0.5 )
		fa = ( indgen( nfreq+1 ) > 1 )^power
		fa(0) = 0
		nf = N_elements( fa )
		if keyword_set( std ) then begin
			fa = ( fa + randomn( sn, nf ) * std ) > 0
			fa(0) = 0
		   endif
		if N_elements( phase ) NE nf then phase = randomu( su, nf )
		fc = fa * exp( 2*!PI * complex(0,1) * phase )
		ft = [ fc(0:nfreq-1), rotate( conj( fc(1:*) ), 2 ) ]
		return, float( fft( ft, 1 ) )
	     END

	  3: BEGIN
		message,"fractal dimension = " + $
			strtrim( ((edim-hpar)>2)<3, 2 ),/INFO
		nf = nfreq + 1
		fa = total( gridxy( nf(0), nf(1), POWER=2 ), 3 )	;squared
		fa(0,0) = 1
		fa = fa^( -0.5 * ( Hpar + 1 ) )		;square-root included.
		fa(0,0) = 0
		fa = [ rotate( fa(1:*,*), 5 ), fa(0:nfreq(0)-1,*) ]
		sz = size( fa )
		if keyword_set( std ) then begin
			fa = ( fa + randomn( sn, sz(1), sz(2) ) * std ) > 0
			fa(sz(1)/2,0) = 0
		   endif
		if N_elements( phase ) NE N_elements( fa ) then $
				phase = randomu( su, sz(1), sz(2) )
		fc = fa * exp( 2*!PI * complex(0,1) * phase )
		return, rfftinv( fc, ft, /CENTER, /FORCE_REAL )
	     END

	  else: BEGIN
		message,"case not implemented",/INFO
		return,(0)
		END
	 ENDCASE
end
