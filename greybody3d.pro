;+
; NAME:
;	greybody3d.pro
; PURPOSE:
;       interate greybody.pro over an [M x N] image cube or sub-cube &
;        stuff output parameters into appropriately-sized arrays
; CALLING:
;	gb = greybody3d, parcube, puncube, corners, modelno...
; INPUTS:
;	parcube = 3D array to contain fit params @ each pixel of image
;	puncube = 3D array of errors in fit params @ each pixel of
;	          image
;       corners = 5-tuple of min & max x & y coords to iterate over, & an
;		origin (0/1) to specify whether the first index is 0 or 1;
;		accounting for the origin is taken care of in iter_gb3d.pro
; COMMON BLOCKS:
;	keys: keylist, par_info ;;from gbpstruct
;	data: nudat, fcube, ecube, wcube;; from iter_gb3d
;	fitcon: Nwav, dof, imd, fwhm, bounds;; also from iter_gb3d
; EXTERNAL CALLS:
;	gbpstruct.pro
;	iter_gb3d.pro <-- see this for more details
;	*Assumes !csi is defined, is in SI units, & contains mass of hydrogen atom
; HISTORY:
;	Written: Rebecca Pitts, 2017.
;---------------------------------------------------

pro greybody3d, parcube, puncube, modcube, corners, modelno

  common keys, keylist, par_info ;;from gbpstruct
  common data, nudat, fcube, ecube, wcube 
  common fitcon, Nwav, dof, imd, fwhm, bounds, kclrcor
  ; note that fwhm = bmd in iter_gb3D
  ; corners has the origin=1/0 discrepancy rectified in iter_gb3d now

  IF corners EQ !null THEN BEGIN
     x0 = 0
     x1 = imd[0]-1
     y0 = 0
     y1 = imd[1]-1
  ENDIF ELSE BEGIN
     x0 = corners[0]
     x1 = corners[2]
     y0 = corners[1]
     y1 = corners[3]
  ENDELSE

  FOR I=x0,x1 DO BEGIN
     FOR J=y0,y1 DO BEGIN
        spec = reform(fcube[I,J,*])
	IF corners EQ !null THEN BEGIN
	   I0 = I
	   J0 = J
	ENDIF ELSE BEGIN
           I0 = I-x0 ;parcube & puncube may be sub-regions
           J0 = J-y0
	ENDELSE

        IF (ecube NE !null) THEN BEGIN
           IF (N_elements(ecube) LT Nwav) THEN err = reform(ecube) ELSE $
              err = reform(ecube[I,J,*])
        ENDIF ELSE err = !null
        IF (wcube NE !null) THEN BEGIN
           IF (N_elements(wcube) LT Nwav) THEN wt = reform(wcube) ELSE $
              wt = reform(wcube[I,J,*])
        ENDIF ELSE wt = !null

        IF (err NE !null) THEN goodi = where(spec GT err, /null) ELSE goodi = where(spec GT bounds, /null)
        IF (modelno NE -1) THEN BEGIN
	   IF (err NE !null) THEN fatali = where(spec[1:-2] LE err[1:-2], count, /null) ELSE $
              fatali = where(spec[1:-2] LE bounds, count, /null)
	ENDIF ELSE count = 0

        IF (N_elements(goodi) LT dof) OR (Nwav - count LT dof) THEN BEGIN ;;count is from fatali
           parcube[I0,J0,*] = !values.d_nan
           puncube[I0,J0,*] = !values.d_nan
           modcube[I0,J0,*] = !values.d_nan

;        ENDIF ELSE IF modelno EQ 1 THEN BEGIN
;	   IF (N_elements(goodi) LT dof+1) OR (Nwav - count LT dof+1) THEN BEGIN ;;count is from fatali
;              parcube[I0,J0,*] = !values.d_nan
;              puncube[I0,J0,*] = !values.d_nan
;              modcube[I0,J0,*] = !values.d_nan
;	   ENDIF

        ENDIF ELSE IF (count GT 0) THEN BEGIN
		;;alternative condition: (spec[0] LT bounds) OR (spec[-1] LT bounds) THEN BEGIN
	   print,"running greybody..."
           IF (N_elements(err) LE 1) THEN errarr = !null ELSE errarr = err[goodi]
           IF (N_elements(wt) LE 1) THEN wtarr = !null ELSE wtarr = wt[goodi]
           gbstruct = greybody(nudat[goodi], spec[goodi], ERRARR=errarr, PAR_INFO=par_info, MODELNO=modelno, $
                               THRESHOLD=bounds, WEIGHT=wtarr, KCLRCOR=kclrcor, /autocall)
           parcube[I0,J0,*] = reform(gbstruct.params)
           puncube[I0,J0,*] = reform(gbstruct.perrs)
	   modcube[I0,J0,goodi] = reform(gbstruct.fmodel)

        ENDIF ELSE BEGIN
           gbstruct = greybody(nudat, spec, ERRARR=err, PAR_INFO=par_info, MODELNO=modelno, $
                               THRESHOLD=bounds, WEIGHT=wt, KCLRCOR=kclrcor, /autocall)
           parcube[I0,J0,*] = reform(gbstruct.params)
           puncube[I0,J0,*] = reform(gbstruct.perrs)
	   modcube[I0,J0,*] = reform(gbstruct.fmodel)

        ENDELSE

        IF (puncube[I0,J0,0] GT 0.5*parcube[I0,J0,0]) OR (puncube[I0,J0,3] GT 0.5d) THEN BEGIN
	   parcube[I0,J0,*] = !values.d_nan
           puncube[I0,J0,*] = !values.d_nan
           modcube[I0,J0,*] = !values.d_nan
	ENDIF
     ENDFOR
  ENDFOR
END

