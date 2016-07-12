;+
; NAME:
;	Mosaic_Main
;
; PURPOSE:
;	Main program for calling:	Prep
;					Mosaic
;					Display_Mosaics
;					Analyze_Mosaics
; EXECUTION:
;	IDL> .run  mosaic_main
;
;	Further execution is all menu-driven.
;	Information, status, and directions are printed in the terminal window
;	and the info/menu_window, and the windows displaying images.
;
; MAIN VARIABLES:
;
;	task = string, next task requested by user.
;
;	menu_window = window # in which informational messages are displayed,
;		holdover from SunView development, since in SunView
;		all menus had to be associated with a window.
;
;	inventory = array of structures containing images, locations, info, etc.
;		created in pro Prep by reading data files.
;
;	SS_FF = array of structures resulting from preprocessing (in Prep),
;		either sky-subtractions (SS) or flat-fielding (FF),
;		can be exported to Mosaic, or saved to a file.
;
;	Raw_Mosaic = array of structures containing images, locations, info, ...
;		imported from Prep to Mosaic, or created in Mosaic.
;
;	mosaic = single image, result of averaging or splicing a Raw_Mosaic.
;	mosaic_spec = structure containing display specifications for mosaic.
;	mosaic_info = the Raw_Mosaic structure minus the images.
;
;	math_result = 2D array (image) result of computations.
;	math_spec = structure defining the math_result.
;
;   common mosaic_array:
;
;	mosaics = array of pooled image data, pointed to by mosaic_specs.
;	mosaic_specs = array of structures containing specs. and
;			pointers into the pool-array: mosaics.
;	mosaic_infos = array of structures containing info about
;			how each mosaic image was formed.
;
;   common math_mosaics:
;
;	mathmos_List = array of structures containing image locations and
;		pointers into the pool-arrays: math_images & math_imscaled.
;	math_images = pool-array of image data, pointed to by mathmos_List.
;	math_imscaled = pool-array image data scaled into bytes for display.
;
;	Elements of math_images are extracted from array mosaics.
;
;   common contour_marks:
;
;	sources = array of structures containing user defined marks to be drawn
;		on images displayed in the "contour mosaic" function,
;		found in the "DISPLAY mosaics" module, or in "ANALYZE mosaics".
;		
;	Labels = array of structures containing user defined Labels
;		to be drawn on images in the "contour mosaic" function.
;
; HISTORY:
;	Written: Frank Varosi, NASA/GSFC, 1990.
;-

if N_elements( select ) NE 1 then begin
	select = 4
   print,"+-------------------------------------------------------------------+"
   print,"+  MOSAIC creation, display and analysis of IR array camera images  +"
   print,"+              by Frank Varosi and Dan Gezari, NASA/GSFC, 1989-2000 +"
   print,"+  ver. 4.4.2, and Frank Varosi at University of Florida, 2001-2012 +"
   print,"+-------------------------------------------------------------------+"
	wait,1
   endif

common mosaic_array, mosaics, mosaic_specs, mosaic_infos
common contour_marks, sources, Labels
common math_mosaics, mathmos_List, math_images, math_imscaled

common control, cur_pro, win_mos, win_prep
if N_elements( win_mos ) NE 1 then win_mos = 4
if N_elements( win_prep ) NE 1 then win_prep = 8

common menus, menu_window

menu_window, XPOS=0, YPOS=!DEVY-128, COL_TAB=col_table, $
             TIT="Info and Instructions for MOSAIC version: 4.4.2"

main_menu = [	"Choose Function:"		,$
		" "				,$
		"PreProcess  images"		,$
		" "				,$
		"MOSAIC creation"		,$
		" "				,$
		"DISPLAY mosaics"		,$
		" "				,$
		"ANALYZE mosaics"		,$
		" "				,$
		"IDL>"				]
task = ""

while (task NE "IDL>") do begin

	CASE task OF

	  "PreProcess":	prep, inventory, SS_FF, Raw_Mosaic, task

	  "MOSAIC":	mosaic, Raw_Mosaic, inventory, SS_FF, mosaic, task

	  "DISPLAY":	display_mosaics, mosaic, mosaic_spec, mosaic_info, task

	  "ANALYZE":	analyze_mosaics, math_result, math_spec, task, deconv_info

	  else: BEGIN
			window_set_show, menu_window, DELAY=0.1
			select = wmenux( main_menu, INIT=select, TITLE=0 )
			task = next_word( main_menu[select>1] )
		  END

       ENDCASE
endwhile

print,"*"
print," to resume MOSAIC type:  IDL> .go"
print,"      or if that fails:  IDL> .run mosaic_main"
task = ""
end
