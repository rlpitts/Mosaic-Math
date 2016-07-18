;+
; NAME:
;	interp_Tau_Eff
; PURPOSE:
;	Compute effective optical depth of a sphere including scattering.
; EXPLANATION:
;	Effective optical depth is estimated by interpolating a grid
;	of Monte Carlo simulations stored in a structure in common.
;	Inputs can be scalar or arrays, but must be all same dimensions.
; CALLING:
;	Tau_Eff_scat = interp_Tau_Eff( gcos, albedos, taus )
; INPUTS:
;	gcos = average < cos( scattering angle ) >
;	albedos = scattering albedos.
;	taus = optical depths (absorbtion + scattering cross-sections).
; OUTPUTS:
;	Function returns scalar/array of effective optical depths
; EXTERNAL CALLS:
;	function interleave
; COMMON BLOCKS:
;	common interp_Tau_Eff, mcrt_gat    ;Monte Carlo rad. trans. simulations.
; PROCEDURE:
;	Interleave the inputs into the Monte Carlo grids and then
;	use the IDL intrinsic function interpolate.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1996.
;-

function interp_Tau_Eff, gcos, albedos, taus

 common interp_Tau_Eff, mcrt_gat

	gp = mcrt_gat(*,0,0).g
	alb = reform( mcrt_gat(0,*,0).albedo )
	taut = reform( mcrt_gat(0,0,*).tau_tot )

	ig = ( interleave( gp, gcos ) > 0 ) < (N_elements( gp )-2)
	ia = ( interleave( alb, albedos ) > 0 ) < (N_elements( alb )-2)
	it = ( interleave( taut, taus ) > 0 ) < (N_elements( taut )-2)

	ng = N_elements( gcos )
	na = N_elements( albedos )
	nt = N_elements( taus )
	np = max( [nt,na,ng] )
	if (ng LT np) then ig = replicate( ig(0), np )
	if (na LT np) then ia = replicate( ia(0), np )
	if (nt LT np) then it = replicate( it(0), np )

	fg = ( gcos - gp(ig) ) / ( gp(ig+1) - gp(ig) ) + ig
	ft = ( taus - taut(it) ) / ( taut(it+1) - taut(it) ) + it
	fa = ( albedos - alb(ia) ) / ( alb(ia+1) - alb(ia) ) + ia

return, interpolate( mcrt_gat.tau_eff, fg, fa, ft )
end
