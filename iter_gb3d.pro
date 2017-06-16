;+
; NAME:
;       iter_gb3d.pro
; PURPOSE:
;	Compute and fit greybody emission spectrum to given flux data
;	at for an array of wavelengths.
;	Note: lots of helper functions ahead of the main one!
;	Uses the following relationships:
;		kappa = kappa_0*(wl_0/wl)^beta
;		kappa_0*N_H = tau_0, where N = col. density & tau << 1
;	If tau = absorption cross-section & N_H = # of H atoms, then
;		kappa*mu*m_H/G2D = tau/N_H = (tau_0/N_H)*(wl_0/wl)^beta
;	     -->tau = N_H*kappa_0*(mu*m_H/G2D)*(nu/nu_0)^beta
;	but that may make N_H unmanageably big, & conversion to those units hard
; CALLING:
;	gb = greybody3D(x, fdat, ... see kwarg list)
;	***Some kwargs can be a float, or 2-, or 3-tuple:
;	arg = init_value -OR-
;	arg = [init_value, bool_fixed] -OR-
;	arg = [init_value, min, max]
;	      (use -1 for min/max to sub in default min/max respectively)
; INPUTS:
;	x = array of wavelengths or frequencies, in any units allowed for XUNIT.
;		defaults to wavelength in microns
;	fluxarr = 3D array of fluxes, i.e. stack of images, at each wavelength
;		in either 'Jy/asec^2', 'MJy/sr', 'Jy/beam', or 'W/m^2/Hz'
; KEYWORDS:
;	ERRCUBE = 1D or 3D array of errors in flux, same units as fluxarr;
;		if 1D, must be same length as x
;		overrides weight if both are present
;	WEIGHT = float array of weights for each flux, nominally in 
;		units of 1/flux; nullified if errcube present
;		see greybody.pro for default
;	XUNIT = string units of abscissa
;		accepts 'm', 'mm', 'um', 'nm', 'Hz', 'kHz', 'MHz', 'GHz', 'THz' 
;		default is 'um'
;	FUNIT = string units of fdat
;		accepts 'MJy/sr','Jy/beam','Jy/asec2', 'W/m2/Hz', or all of the
;		above without the '/'
;	THRESHOLD = lower limit of valid flux values; default is 0.
;		may incorporate upper limits later
;	BGBOX = can be a 4-tuple, [x0,y0,x1,y1] defining bounds of area
;	        considered background to be subtracted, or a string to pick 
;               out the row containing those coords in BGSubRegions.tbl (a 
;               look-up table of background regions), or a 2-tuple in
;               [lookupfile.ext,string] format where .ext can be tbl, dat,
;               or txt, and string is a key to find the row w/in lookupfile.
;               Default lookup table is BGSubRegions.tbl.
;               If not set, BG-subtraction will be skipped entirely.
;       MAPP = Boole; if true, call pmap_maker & make parameter maps
;	SAVEM = Boole; if set, save parameter maps
;       SAVEF = Boole; if set, save params & their errors to fits files
;       SAVEP = Boole; if set, save variables & common blocks to .sav file
;		(warning - these are huge & eat through your quota like candy)
;       SAVEDIR = string path to directory where you want to save
;               data/images; default is '~/mosaic/scratch/'
;       TAIL = string to prepend to map file names so you recognize
;	        them; program auto. specifies parameter being mapped
;	NORMA = kludge factor array: if imcomb fails to normalize pixels
;		by number of input images, divide images depthwise by this
; OUTPUTS:
;	RESULTS: a struct containing cubes of fitting params & their errors @ each pixel
;	Also can save parameter maps & error maps - still working on 
; COMMON BLOCKS:
;	None b/c I'm a Python native & these things bug the scheisse out of me
;	(I kid - it's a neat trick but I'm writing most of these
;	modules to work alone OR together)
; EXTERNAL CALLS:
;	nuFnu2SI.pro
;	gbpstruct.pro
;	greybody.pro & its helper functions
;	*Assumes !const is defined, is in SI units, & contains mass of hydrogen atom
; HISTORY:
;	Written: Rebecca Pitts, 2016.
;	Jan 2017 - updated to run with nuFnu2SI.pro
;		helper fxns x2nu and fucon decomissioned
;	Feb 2017 - added error cube & weight options, plotmap routines
;
;---------------------------------------------------

pro greybody3D, parcube, puncube, corners, modelno

  common keys, keylist, par_info ;;from gbpstruct
  common data, nudat, fcube, ecube, wcube
  common fitcon, Nwav, dof, imd, fwhm, bounds
  ; note that fwhm = bmd in iter_gb3D

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
        I0 = I-x0-1 ;parcube & puncube may be sub-regions
        J0 = J-y0-1

        IF (ecube NE !null) THEN BEGIN
           IF (N_elements(ecube) LE Nwav) THEN err = reform(ecube) ELSE $
              err = reform(ecube[I,J,*])

        ENDIF ELSE err = !null
        IF (wcube NE !null) THEN BEGIN
           IF (N_elements(wcube) LE Nwav) THEN wt = reform(wcube) ELSE $
              wt = reform(wcube[I,J,*])
        ENDIF ELSE wt = !null

        IF (err NE !null) THEN goodi = where(spec GT err, /null) ELSE goodi = where(spec GT bounds, /null)
        IF (err NE !null) THEN fatali = where(spec[1:-2] LE err[1:-2], count, /null) ELSE $
           fatali = where(spec[1:-2] LE bounds, count, /null)

        IF (N_elements(goodi) LE dof) OR (Nwav - count LE dof) THEN BEGIN ;;count is from fatali
           parcube[I0,J0,*] = !values.d_nan
           puncube[I0,J0,*] = !values.d_nan

        ENDIF ELSE IF (count GT 0) THEN BEGIN
		;;alternative condition: (spec[0] LT bounds) OR (spec[-1] LT bounds) THEN BEGIN
           IF (N_elements(err) LE 1) THEN errarr = !null ELSE errarr = err[goodi]
           IF (N_elements(wt) LE 1) THEN wtarr = !null ELSE wtarr = wt[goodi]
           print, N_elements(nudat[goodi]), N_elements(fcube[I,J,goodi]), N_elements(errarr[goodi])
           gbstruct = greybody(nudat[goodi], spec[goodi], ERRARR=errarr, PAR_INFO=par_info, $
                               MODELNO=modelno, THRESHOLD=bounds, WEIGHT=wtarr, /autocall)
           parcube[I0,J0,*] = reform(gbstruct.params)
           puncube[I0,J0,*] = reform(gbstruct.parerrs)

        ENDIF ELSE BEGIN
           gbstruct = greybody(nudat, spec, ERRARR=err, PAR_INFO=par_info,MODELNO=modelno, $
                               THRESHOLD=bounds, WEIGHT=wt, /autocall)
           parcube[I0,J0,*] = reform(gbstruct.params)
           puncube[I0,J0,*] = reform(gbstruct.parerrs)

        ENDELSE
     ENDFOR
  ENDFOR
END

function iter_gb3d, x, datacube, ERRCUBE=errcube, WEIGHT=weight, XUNIT=xunit, FUNIT=funit, $
                    FWHM=fwhm, THRESHOLD=threshold, MODELNO=modelno, BGBOX=BGBOX, MAPP=mapp, $
                    SAVEDIR=savedir, TAIL=tail, HDR=hdr, SAVEF=savef, SAVEP=savep
		    ;;, PAR_INFO=par_info, VERBOSE=verbose

  common keys, keylist, par_info ;;from gbpstruct.pro
  ;;rest assigned piecemeal - some renaming gymnastics required
  common data, nudat, fcube, ecube, wcube 
  common fitcon, Nwav, dof, imd, bmd, bounds
  ;common savecon, root, leaf

  Nwav = N_elements(x)
  dof = N_elements(where(par_info[*].fixed EQ 0, /null))
  IF (dof GT Nwav) THEN MESSAGE, "Error: insufficient data to constrain this many free parameters"
  ;woodstock doesn't have orderedhash capabilities >:(

  IF (isa(datacube,'string') EQ 1) THEN imcube = readfits(datacube,hdr) ELSE $
  IF (isa(datacube,/array) EQ 1) THEN imcube = datacube

  IF (Nwav LE 0) OR (N_elements(imcube) LE 0) THEN RETALL, 'FatalError: Missing required input (check documentation)'
  ;;/dimensions gives size along each axis, /n_dimensions counts total dimensions
  imd = size(imcube,/dimensions) ;;imd[-1] will be the same whether 1D or 3D - I checked
  IF (Nwav NE imd[-1]) THEN BEGIN
     MESSAGE,"Error: unequal image cube depth and length of wavelength array" + string(7b),/INFO
     RETALL,Nwav,imd
  ENDIF

  IF n_elements(xunit) LT 1 THEN xunit='um'
  IF n_elements(funit) LT 1 THEN funit='Jy/beam'

  IF (N_elements(fwhm) LE 0) AND (where(strmatch(['Jy/beam','Jybeam'],funit,/fold_case)) NE -1) THEN BEGIN
     IF (N_elements(hdr) NE 0) THEN BEGIN
	IF (sxpar(hdr,'bmaj') EQ sxpar(hdr,'bmin')) THEN fwhm = sxpar(hdr,'bmaj')*3600.d $ 
	  ELSE fwhm = [sxpar(hdr,'bmaj')*3600.d,sxpar(hdr,'bmin')*3600.d]
     ENDIF ELSE fwhm = 37.d ;;Mopra resolution - won't hurt if funit isn't in Jy/beam
  ENDIF
  IF (where(strmatch(['Jy/beam','Jybeam'],funit,/fold_case)) EQ -1) THEN fwhm = !NULL
  bmd = fwhm ;renaming gymnastics

  IF (strpos(funit,'pixel') NE -1) THEN BEGIN
     IF (N_elements(hdr) NE 0) THEN BEGIN
	CASE 1 OF
	   (sxpar(hdr,'cdelt1') NE 0): pxs = [sxpar(hdr,'cdelt1')*3600.d, sxpar(hdr,'cdelt2')*3600.d]
	   (sxpar(hdr,'cd1_1') NE 0): pxs = [sxpar(hdr,'cd1_1')*3600.d, sxpar(hdr,'cd2_2')*3600.d]
	   (sxpar(hdr,'PXSCALE1') NE 0): pxs = [sxpar(hdr,'pxscale'), sxpar(hdr,'pxscale')]
	   (sxpar(hdr,'PLATESCA') NE 0): pxs = sxpar(hdr,'platesca')
	ENDCASE
     ENDIF ELSE pxs = 12.d ;;Mopra pixel scale - won't hurt if funit isn't in Jy/pixel
  ENDIF
  IF (strpos(funit,'pixel') EQ -1) THEN pxs = !NULL

  SIdata = nufnu2si(x,imcube,xunit,funit,beam=fwhm,pxsz=pxs)
  nudat = SIdata.nu
  fcube = SIdata.flux

  IF (N_elements(errcube) GT 0) THEN BEGIN

     IF (isa(errcube,'string') EQ 1) THEN errcube = readfits(errcube)
     erd = size(errcube,/dimensions)
     IF (total(erd NE imd) NE 0) THEN BEGIN
        MESSAGE,"Error: image cube and error cube dimensions do not match" + string(7b),/INFO
        RETALL,imd,erd
     ENDIF

     IF (N_elements(errcube) EQ Nwav) THEN BEGIN
	errcube = errcube/norma
	dummy = make_array(imd)
	FOR i=0,imd[-1]-1 DO dummy[*,*,i] = errcube[i]
	errcube = dummy
     ENDIF

     ecube = nuFnu2SI(x,errcube,xunit,funit,beam=fwhm,pxsz=pxs,/xoff)
     IF (N_elements(weight) GT 0) THEN BEGIN
        PRINT, 'Error array takes precedence over weight array; weight nullified'
        wcube=!null
     ENDIF

  ENDIF ELSE IF (N_elements(errcube) EQ 0) AND (N_elements(weight) GT 0) THEN BEGIN

     IF (isa(weight,'string') EQ 1) THEN weight = readfits(weight)
     wrd = size(weight,/dimensions)

     IF (wrd[-1] NE imd[-1]) THEN BEGIN
        MESSAGE,"Error: image cube and weight cube dimensions do not match" + string(7b),/INFO
        RETALL,imd,wrd
     ENDIF

     IF (N_elements(weight) EQ Nwav) THEN BEGIN
	dummy = make_array(imd)
	FOR i=0,imd[-1]-1 DO dummy[*,*,i] = weight[i]
	wcube = dummy
     ENDIF ELSE wcube = weight ;common block naming gymnastics
     ecube = !null

  ENDIF ELSE IF (N_elements(errcube) EQ 0) AND (N_elements(weight) EQ 0) THEN BEGIN
     wcube = !null
     ecube = !null
     PRINT, 'Warning: defaulting to internal peak-biased weight scheme.'
     PRINT, 'Do not trust chi-squared values or parameter errors, if they print.'
  ENDIF

  IF n_elements(threshold) EQ 0 THEN threshold = 0.D
  IF (threshold NE 0) THEN bounds = nuFnu2SI(x,threshold,xunit,funit,beam=fwhm,pxsz=pxs,/xoff) ELSE bounds = 0.D

  IF (n_elements(modelno) EQ 0) THEN BEGIN
     IF (N_elements(par_info) LE 6) THEN modelno = 1
     IF (N_elements(par_info) EQ 7) THEN modelno = 2
     IF (N_elements(par_info) EQ 8) THEN modelno = 3
     ;;require explicit selection of modelno = 3 for 
     IF (N_elements(par_info) GT 8) THEN MESSAGE, "Error: too many parameters to fit"
  ENDIF
  
  IF (N_elements(bgbox) EQ 0) THEN BEGIN 
     verts = !null
  ENDIF ELSE BEGIN

     IF (N_elements(bgbox) EQ 1) and (isa(bgbox, 'string') EQ 1) THEN BEGIN

        ;either I'm terribly spoilt or this is a rly dumb way to read files
        OPENR, lun, '~/idl/vlibm/rtfit/BGSubRegions.tbl', /GET_LUN
        ; Read one line at a time, saving the result into array
        array = ''
        line = ''
        WHILE NOT EOF(lun) DO BEGIN
           READF, lun, line
           array = [array, line]
        ENDWHILE
        ; Close the file and free the file unit
        FREE_LUN, lun
        row = array[where(strmatch(array,bgbox+'*',/fold_case))]
        rowsplit = strsplit(row,/extract)
        verts = float(rowsplit[1:4]) ;corners, [x0,y0,x1,y1]

     ENDIF ELSE IF (N_elements(bgbox) EQ 2) THEN BEGIN

        IF (file_search(bgbox[0],/test_read) NE '') THEN BEGIN 
           OPENR, lun, bgbox[0], /GET_LUN
           array = ''
           line = ''
           WHILE NOT EOF(lun) DO BEGIN
              READF, lun, line
              array = [array, line]
           ENDWHILE
        ; Close the file and free the file unit
           FREE_LUN, lun
           row = array[where(strmatch(array,bgbox[1]+'*',/fold_case))]
           rowsplit = strsplit(row,/extract)
           verts = float(rowsplit[1:4])     ;corners, [x0,y0,x1,y1]
        ENDIF ELSE MESSAGE, 'Error: first entry of 2-tuple BGBOX must be a file name.'

     ENDIF ELSE IF (N_elements(bgbox) EQ 4) THEN verts = float(bgbox) ELSE MESSAGE, 'Error: check form of BGBOX'

     FOR k=0,imd[-1]-1 DO BEGIN
	subf = fcube[ verts[0]:verts[2], verts[1]:verts[3], k]
	;; don't subtract unless there's something reliable to subtract
	IF N_elements(where(subf NE !values.d_nan)) NE 1 THEN BEGIN
	   IF min(subf) GT 0 then m = min(subf) ELSE m = median(subf)
           stdv = stddev(subf, /nan)
	   IF (ecube NE !null) THEN ecube[*,*,k] = ecube[*,*,k] + stdv
	   IF (wcube NE !null) THEN wcube[*,*,k] = wcube[*,*,k] + 1/stdv
	   IF (m LE stdv) THEN fcube[*,*,k] = fcube[*,*,k] - m ;; ELSE fcube[*,*,k] = fcube[*,*,k] - stdv
	ENDIF
     ENDFOR
     fcube[where(fcube LT 0)] = 0.0
     
     verts = !null
     ; roughly approximate increase in error due to BG subtraction
     IF (ecube NE !null) THEN ecube = ecube*sqrt(2)
     IF (wcube NE !null) THEN wcube = wcube*sqrt(2)

  ENDELSE

  parcube = fltarr(imd[0],imd[1],N_elements(par_info))
  puncube = fltarr(imd[0],imd[1],N_elements(par_info))
  greybody3D, parcube, puncube, verts, modelno ;by here verts=!null
  IF (total(parcube EQ par_info.value[*]) EQ N_elements(parcube)) THEN $
     MESSAGE,'ValueError: mpfitfun failed over whole image cube'
  IF (N_elements(ecube) EQ 0) THEN results = {nu:nudat,flux:fcube,params:parcube,perrs:puncube} $
     ELSE results = {nu:nudat,flux:fcube,errflux:ecube,params:parcube,perrs:puncube}
  
  IF (N_elements(savedir) EQ 0) THEN savedir='~/mosaic/save'
  IF (N_elements(tail) EQ 0) THEN BEGIN
     pfixnum = 'pfix'+strjoin(strcompress(string(transpose(par_info.fixed[*])),/remove_all))
     IF (isa(bgbox, 'string') EQ 1) THEN tail = bgbox+'_'+pfixnum ELSE tail = pfixnum
  ENDIF
  IF ( strpos(savedir,'/',/reverse_search) EQ strlen(savedir)-1 ) THEN fname = savedir+tail ELSE $
     fname = savedir+'/'+tail
  IF keyword_set(savep) THEN save,/comm,/variables,filename=fname+".sav"

  ;last of common block renaming gymnastics
  ;root = savedir & leaf = tail

  IF keyword_set(mapp) or keyword_set(savef) THEN BEGIN
    ;;IF keyword_set(savef) THEN BEGIN
      pmap_maker, results, HDR=hdr, SAVEFITS=savef, SAVEPLOT=saveplot, SAVEDIR=savedir, TAIL=tail
    ;;ENDIF ELSE pmap_maker, results, HDR=hdr, SAVEDIR=savedir, TAIL=tail
  ENDIF

  RETURN, results
END


