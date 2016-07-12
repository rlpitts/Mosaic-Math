;+
; FILE:
;	mosaic_setup.pro
; PURPOSE:
;	Do NOT modify this file!
;	This file should be linked to: ~/mosaic/setup.pro
;	so that it is automatically executed by ~/idl/idl_startup.pro
;	(need:  setenv  IDL_STARTUP  ~/idl/idl_startup.pro).
;	Copy and put any desired modifications and variable settings in: 
;	~/mosaic/setup_user.pro,  which is invoked at end of this IDL script.
; HISTORY:
;	Written by Frank Varosi, NASA/GSFC 1989.
;-

!PATH = "~/mosaic/code/lib:"		+ $
	"~/mosaic/code/lib/image:"	+ $
	"~/mosaic/code/lib/struct:"	+ $
	"~/mosaic/code/prep:"		+ $
	"~/mosaic/code/mosaic:"		+ $
	"~/mosaic/code/display:"	+ $
	"~/mosaic/code/analyze:"	+ $
	"~/mosaic/code:"		+ !PATH

; The following common block variables are the default setups for
; the MOSAIC image processing software.
; The variables can be redefined in setup_user.pro (then overrides these defs.)
; Most directories are specified as string arrays giving path in dir-tree,
; independent of op.sys., the file system syntax is applied later in code.
; For example: ["data","raw_mosaics"] becomes "data/raw_mosaics" in UNIX.

common dir_names, dirinv, dirmos_raw, dirmos_aver	;sub-directories for 
dirinv = ["data","inventory"]				; storing results.
dirmos_raw = ["data","raw_mosaics"]
dirmos_aver = ["data","averaged_mosaics"]

common dir_name1, dir_FITS
dir_FITS = ["data","fits"]	;for Load FITS option in display_mosaics

common dir_name2, dirmos_math
dirmos_math = ["data","math_mosaics"]	;for pro analyze_mosaics

common select_files, dir_obs		;directory of observation data files.
dir_obs = "/data/obs/"			;must terminate with / for unix.

common array_scale, ArcSec_Pix_x, ArcSec_Pix_y	   ;plate scale of pixel array
ArcSec_Pix_x = .26				;can also be set in menu option:
ArcSec_Pix_y = .26				; "mosaic by offsets",
				;or set in contour with "plate scale" option.

common color_options, ct_reset, ct_feedback
common color_table0, col_table

col_table = 3	;initial IDL color table #.
ct_reset = 0	; 0 = do not reset adjct transfer function after loading.

common mosaic_options, border_images, border_groups

border_images = 0	;border with black on pop (menu: set mosaic defaults)
border_groups = 255	;border with max (white) on pop group.

common display_options, Magf_min, Magf_max
common display_option2, rotation_default

Magf_min = 0		;set to 0 if Magf < 1 is desired (in setup_user.pro).
rotation_default = 1	;set to 0 for no rotation.

common contour_init, site_name	;set site_name="your place" in setup_user.pro

device,RETAIN=2

@setup_user		;invoke site specific override of setups...
