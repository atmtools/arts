; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
PRO showme, jobname=jobname,          $
            jobdir=jobdir,            $
            plotfilename=plotfilename
;
;==========================================================================
;
; NAME:
;       showme
; PURPOSE:
;       provide a simple call routine to plot the absorption per tag group
; EXPLANATION:
;       This IDL procedure just calls the procedure 
;
; CALLING EXAMPLES:
;       showme, jobname='myartscontrolfilename',  $
;               jobdir ='~/arts/',                $
;               plotfilename='myabsplot'
;
; INPUTS:
;       jobname      (string)  name of the arts job from which the
;                              output should be displayed
;       jobdir       (string)  directory wher the output of the arts
;                              job is located
;       plotfilename (string)  name of the postscript file
;
; OUTPUTS:
;       Postscript file with absorption plot.
;
; CALLS:
;       plot_abs_per_tg
;
; MODIFICATION HISTORY:
;       04/27/01  TKS  alpha version created 
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
; input control:
;
if not keyword_set(jobname) then begin
    jobname = 'cont'
endif
;
if not keyword_set(jobdir) then begin
    spawn,'pwd', jobdir
    jobdir=jobdir+'/'
endif
;
if not keyword_set(plotfilename) then begin
    plotfilename = 'cont_plot_dBkm'
endif
;
print, 'arts jobname   : ',jobname
print, 'arts jobdir    : ',jobdir
print, 'output filename: ',plotfilename
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
plot_abs_per_tg, jobname,                    $
                 f,                          $
                 abs,                        $
                 alt,                        $
                 altitude=10.0,              $
                 jobdir=jobdir,              $
                 pressure='hPa',             $
                 temperature='K',            $
                 absunit='dB/km',            $
                 plotfilename=plotfilename,  $
                 plotfileformat=2,           $
                 plotsum=1
;
ende:
END
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
