; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
PRO showme, jobname=jobname,          $
            jobdir=jobdir,            $
            plotfilename=plotfilename,$
            altitude=altitude
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
;       showme, jobname='wats_pr1', jobdir='/home/home01/tkuhn/ARTS/'
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
    plotfilename = jobname
endif
;
if not keyword_set(altitude) then begin
    selalt = 1.0
endif else begin
    selalt = altitude
endelse
print, '----------------------------------------------------------'
print, 'arts jobname     : ',jobname
print, 'arts jobdir      : ',jobdir
print, 'output filename  : ',plotfilename
;
spawn, 'rm -f '+plotfilename+'.ps'
print, ' '
print, 'selected altitude: ',selalt,' km'
print, '----------------------------------------------------------'
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
goto, absplot
;
;
vmrplot:
pvmrname = plotfilename+'_vmr'
plot_vmr_per_tg, jobname=jobname,            $
                 jobdir=jobdir,              $
                 punit='hPa',                $
                 tunit='K',                  $
                 vmrunit='VMR',              $
                 altunit='km',               $
                 plotfilename=pvmrname,      $
                 plotfileformat=1
;goto, ende

absplot:
pabsname = plotfilename+'_abs'
plot_abs_per_tg, jobname,                    $
                 f,                          $
                 abs,                        $
                 alt,                        $
                 altitude=selalt,            $
                 jobdir=jobdir,              $
                 pressure='hPa',             $
                 temperature='K',            $
                 absunit='1/km',             $
                 plotfilename=pabsname,      $
                 plotfileformat=2,           $
                 plotsum=2,           $
                 avoid_tg=["N2"] 
;
ende:
END
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
