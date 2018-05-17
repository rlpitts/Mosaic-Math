function bg_model, inmap, uinmap, MOMAP=momap, ADDMASK=addmask, UNMASK=unmask, $
	MASKFILE=maskfile, OUTMAP=outmap, THRESHOLD=threshold, SAVEBG=savebg

;+
;Name:
;	bg_model
;Purpose:
;	Use source map to mask image cube, bilinearly interpolate over
;	masked values to model off-source emission in every wavelength
;	slice of the image cube, subtracts model from image, & saves
;	result (& optionally background model map)
;Calling Sequence:
;	ofname = bg_model(inmap,MOMAP=momap,ADDMASK=[[x1,x2,...],[y1,y2,...]], $
;				OUTMAP=outmap,THRESHOLD=threshold,/SAVEBG)
;Inputs:
;	INMAP = str name of image cube to model backgrond of & subtract from
;	UINMAP= str name of uncertainty cube of INMAP
;	MOMAP = str name or list of source map name(s)
;	ADDMASK = 2xN array of X & Y coords of corners of polygonal region to
;		mask in addition to the parts of MOMAP above THRESHOLD
;	UNMASK = 2xN array of X & Y coords of corners of polygonal region to
;		UNmask - applied after masking is done
;	MASKFILE = ds9.reg file with coordinates of regions to +/- mask; assumes that
;		1) all regions have distinct names & 2) polygons with a hole are given 
;		by 2 regions with the inner polygon named after the outer with some 
;		keyword containing "in" (e.g. "_inner", "_interior", "_inside", etc)
;Optional Inputs:
;	OUTMAP = name of file to write background-subtracted image cube to;
;		defaults to inmap:r+"_bgsubbed.fits" (to use the csh shorthand);
;		Uncertainty cube will be named inmap:r+"_err.fits"
;	THRESHOLD = lower bound of what constitutes "source" emission in MOMAP;
;		defaults to 10th percentile
;	SAVEBG = keyword to prompt program to save background model cube to file;
;		file will be saved in same directory as inmap, named 
;		inmap:r + "_bgmap.fits" (to use the csh shorthand)
;Outputs:
;	Name of output file, for convenience
;External Calls:
;	wheretomulti, interpolate, sxaddpar, writefits, cgsetdifference
;-

   indat = readfits(inmap,hdr)
   bgdat = readfits(inmap)
   indat = indat+min(indat) & bgdat = bgdat+min(bgdat) ;;(result will be relative)
   uindat = readfits(uinmap)
   ubgdat = readfits(uinmap)
   ;;need copies to subtract from
   nm = n_elements(momap)
   imd = size(indat,/dimensions)
   umd = size(uindat,/dimensions)
   if total(imd) NE total(umd) then begin
      MESSAGE, "Error: image cube and error cube must have the same dimensions" &
      STOP
   endif
    ;;^this should just never happen, but you know what they say about idiot-proofing

   if nm GT 1 then begin
      mofi = readfits(momap[0])
      for j=1,nm-1 do begin
         mofj = readfits(momap[j])
	  ;;don't use mean; should be little if any overlap
	  ;;if there is overlap, it should be summed b/c
	  ;;dust map doesn't distinguish v components
         mofi = total([[[mofi]],[[mofj]]],3,/double,/nan)
      endfor
   endif else mofi = readfits(momap)

   modim = size(mofi,/dimensions)
   if modim[0] NE imd[0] or modim[1] NE imd[1] then begin
      MESSAGE, "Image cube & Mopra map must have the same first 2 dimensions"
      STOP
   endif

   mofi[where(mofi LE 0,/null)] = !values.d_nan
   ms = mofi[sort(mofi)] & msf = ms[where(finite(ms))]
    ;;I set threshold = mean - stdev --> too low, but even 0.75 sigma was too high
    ;; otherwise mask 10th percentile & up (I liked where the ds9 contours landed)
   if n_elements(threshold) LE 0 then threshold=msf[ceil(n_elements(msf)/10)]
   maskn = where((mofi GT threshold))

   if N_elements(addmask) NE 0 then begin
      check = size(addmask)
      if (check[0] NE 2) then begin
	 MESSAGE, "Error: ADDMASK must be in the form [[x1,x2,...],[y1,y2,...]]" & STOP
      endif
      if (check[1] LT 3) then begin
	    MESSAGE, "Error: polygon specified by ADDMASK must have at least 3 vertices" & STOP
      endif
      temp = [ maskn, polyfillv(addmask[*,0],addmask[*,1],imd[0],imd[1]) ]
      maskn = temp[uniq( temp,sort(temp) )]
   endif

   if N_elements(maskfile) NE 0 then begin
      if strmatch(maskfile,'*.reg',/fold_case) then begin
	 OPENR, lun, maskfile, /get_lun
	 strarr = [] & line = ''
	 WHILE NOT EOF(lun) DO BEGIN
	    READF, lun, line
            if strmatch(line,'polygon(*',/fold_case) then strarr = [strarr, line]        
	 ENDWHILE
	 FREE_LUN, lun
	 remsubs=[]
	 foreach row,strarr,i do begin
	    vals = strsplit( strmid( row,strpos(row,'(')+1,strpos(row,')')-strpos(row,'(')-1 ),',',/extract )
	    flag = strmid(line,strpos(line,'{'),strpos(line,'}')-strpos(line,'{'))
	    mx=[floor(double(vals[0:*:2]))]
	    my=[floor(double(vals[1:*:2]))]
	    if total(strmatch(flag,'*in*',/fold_case)) EQ 1 then begin
		remsubs = [remsubs,polyfillv(mx,my,imd[0],imd[1])] ;'inner','inside',etc flag interior annuli bounds
	    endif else begin
		temp = [ maskn, polyfillv(mx,my,imd[0],imd[1]) ]
		maskn = temp[uniq( temp,sort(temp) )]
	    endelse
	 endforeach
	 if n_elements(remsubs) NE 0 then maskn = cgsetdifference(maskn,remsubs)
      endif ;;else begin...I'll add this another time, maybe
   endif

   if N_elements(unmask) NE 0 then begin
      check = size(unmask)
      if (check[0] NE 2) then begin
	 MESSAGE, "Error: UNMASK must be in the form [[x1,x2,...],[y1,y2,...]]" & STOP
      endif
      if (check[1] LT 3) then begin
	    MESSAGE, "Error: polygon specified by UNMASK must have at least 3 vertices" & STOP
      endif
      remsubs = [remsubs,polyfillv(unmask[*,0],unmask[*,1],imd[0],imd[1])]
      maskn = cgsetdifference(maskn,remsubs)
   endif
      
   wheretomulti,mofi,maskn,mcol,mrow
    ;;mcol = maskn mod imd[0]
    ;;mrow = maskn / imd[0]

   for j=0,(n_elements(mcol)-1) do bgdat[mcol[j],mrow[j],*] = !values.d_nan
   for j=0,(n_elements(mcol)-1) do ubgdat[mcol[j],mrow[j],*] = !values.d_nan
   for k=0,imd[2]-1 do begin ;;interp fails, blank area is too big
      slice=bgdat[*,*,k]
      uslice=ubgdat[*,*,k]
      refpts=where(finite(bgdat[*,*,k]))
      wheretomulti,bgdat[*,*,k],refpts,Xi,Yi
      triangulate,Xi,Yi,tris,bpts
      ivals=griddata(Xi,Yi,slice[refpts],/linear,triangles=tris,xout=mcol,yout=mrow)
      check=where(finite(ivals),count)
      ;print,count
      ;print,size(ivals),ivals[1:10]
      uivals=griddata(Xi,Yi,uslice[refpts],/linear,triangles=tris,xout=mcol,yout=mrow)
      for j=0,(n_elements(mcol)-1) do bgdat[mcol[j],mrow[j],k]=ivals[j]
      for j=0,(n_elements(mcol)-1) do ubgdat[mcol[j],mrow[j],k]=uivals[j]
      if (where(~finite(bgdat[*,*,k]),/null) NE !null) then begin
	 MESSAGE,"Error: interpolation referenced invalid data" & STOP
      endif
   endfor
      ;;IDL makes it a royal PITA to insert interpolated values into original array
      ;; b/c it allows (possibly expects) interpolation to non-integer subscripts
   sxaddpar,hdr,'HISTORY','Background model generated with bg_model.pro'
   
   if (where(~finite(bgdat),/null) NE !null) then begin
      MESSAGE,"Error: check variables ibgdat and interplocs" & STOP
   endif else print,'Bilinear interpolation complete.'

   if keyword_set(savebg) then begin
      bhdr = hdr
      sxaddpar,bhdr,'HISTORY','to be subtracted from '+inmap
      s = strsplit(inmap,'.',/extract)
      bgmap = s[0]+'_bgmap.'+s[-1]
      ubgmap = s[0]+'_bgmap_err.'+s[-1]
      writefits, bgmap, bgdat, bhdrpritn
      writefits, ubgmap, ubgdat, bhdr
      sxaddpar,hdr,'HISTORY','background model: '+bgmap
      sxaddpar,hdr,'HISTORY','background model uncertainty: '+ubgmap
   endif
   
   outdat = indat - bgdat
   outdat[where(outdat LT 0)] = 0.d
   if total(outdat,/nan) EQ 0 then begin
      MESSAGE,"Error: outdat = bgdat" & STOP
   endif else print,'Writing results to file...'
   outunc = sqrt( uindat^2 + ubgdat^2 )

   if n_elements(outmap) LE 0 then begin
      s = strsplit(inmap,'.',/extract)
      outmap = s[0]+'_bgsubbed.'+s[-1]
   endif
   s2 = strsplit(outmap,'.',/extract)
   outumap = s2[0]+'_err.'+s2[-1]
   sxaddpar,hdr,'HISTORY','original image cube: '+inmap
   if nm GT 1 then sxaddpar,hdr,'HISTORY','source maps: '+strjoin(momap,'; ') $
   else sxaddpar,hdr,'HISTORY','source map: '+momap
   writefits, outmap, outdat, hdr
   sxaddpar,hdr,'HISTORY','original error cube: '+uinmap
   writefits, outumap, outunc, hdr

   return,[outmap, outumap]
end
