;; ==========================================================================
;; ####################### ARTS IDL INTERFACE PROCEDURE #####################
;; ==========================================================================
;;
;; NAME    : aii_compile
;;
;; PURPOSE : this batch file sets the frame for an ARTS IDL session using
;;           the functions and procedures of the aii directory and some
;;           external libraries 
;;
;; CALLING : in the IDL session type '@aii_compile'
;;
;; EXTERNAL: external libraries are 
;;           TeXtoIDL (http://physweb.mnstate.edu/mcraig/TeXtoIDL/)
;;           astron   (http://idlastro.gsfc.nasa.gov/homepage.html)
;;           idlps    (ftp.astro.washington.edu)
;;           users not from iup-sat need to download the files and to 
;;           upodate the !PATH global variable below
;; 
;; HISTORY : alpha verison 2003-04-03, Thomas Kuhn, iup Bremen
;;
;;
;; #########################################################################

print, ' '
print,'=====================< start of AII_COMPILE >====================='
print, ' '

;; Student t-distribution (statistics)
.COMPILE aii_student_t_dist_table

;; file checks
.COMPILE aii_checks

;; retrieves the host name (UNIX)
.COMPILE aii_get_hostname

;; retrieves the datum from the computer clock
.COMPILE aii_writedatum

;; read/write an ARTS variable
.COMPILE read_artsvar
.COMPILE write_artsvar

;; returns the tag groups found in an ARTS controlfile
.COMPILE read_tag_groups

;; rearrange the absorption and tag groups according to their absorption magnitude
.COMPILE sort_abs_2
.COMPILE sort_abs

;; read/write data from/to a file in ARTS format
.COMPILE write_datafile
.COMPILE read_datafile
.COMPILE aii_aa_read_general

;; plotting functions/procedures
.COMPILE aii_color_table
.COMPILE aii_klegend_d
.COMPILE aii_prologue_l
.COMPILE aii_plotsymbols
.COMPILE aii_epilogue
.COMPILE aii_plot_legend
.COMPILE aii_plot_file

;; plotting absorption, TB, transmission or VMRs
.COMPILE plot_abs_per_tg
.COMPILE plot_abs_per_tg_2
.COMPILE plot_RT_1
.COMPILE plot_trans_1
.COMPILE plot_vmr_per_tg
.COMPILE showme

;; creating an ARTS control file
.COMPILE CreateArtsControlFile


;; fit of H2O continuum parameters
.COMPILE WVContParamFit


;; water vapor saturation pressure calculation
.COMPILE WaterVaporSatPressure
.COMPILE RH_test


print, ' '
print,'======================< end of AII_COMPILE >======================'
print, ' '

;; #########################################################################
