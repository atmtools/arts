; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
PRO showme, jobname=jobname,          $
            jobdir=jobdir,            $
            run=run,                  $
            plotfilename=plotfilename,$
            altitude=altitude,        $
            avoid_tg=avoid_tg 
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
close, /all
;
; input control:
;
if not keyword_set(jobname) then begin
    jobname = 'O2a'
endif
;
if not keyword_set(jobdir) then begin
    spawn,'pwd', jobdir
    jobdir='~/PhD/tksthesis/figures/O2/'
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


if not keyword_set(avoid_tg) then begin
    avoid_tg = ''
endif


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
if not keyword_set(run) then begin
    print,'possible keywords for run:'
    print,'  >   VMR'
    print,'  >   ABS'
    print,'  >   RT'
    print,'  >   TRANS'
    goto, ende
endif
if (STRUPCASE(run) EQ 'VMR')   then goto, vmrplot
if (STRUPCASE(run) EQ 'ABS')   then goto, absplot
if (STRUPCASE(run) EQ 'RT')    then goto, RTplot
if (STRUPCASE(run) EQ 'TRANS') then goto, TRAplot
;
; --------------------------------------------------------------------------
;
vmrplot:
pvmrname = plotfilename+'_vmr'
plot_vmr_per_tg, jobname=jobname,            $
                 jobdir=jobdir,              $
                 punit='hPa',                $
                 tunit='K',                  $
                 yunit='ALTITUDE',           $
                 vmrunit='VMR',              $
                 altunit='km',               $
                 plotxaxis='LOG',            $
                 plotyaxis='LIN',            $
                 yrange=[0.0,15.0],          $
                 plotfilename=plotfilename,  $
                 plotfileformat=1,           $
                 avoid_tg=avoid_tg
goto, ende

; --------------------------------------------------------------------------

absplot:
pabsname = plotfilename+'_abs'
plot_abs_per_tg_2, jobname,                  $
                 f,                          $
                 abs,                        $
                 alt,                        $
                 altitude=selalt,            $
                 jobdir=jobdir,              $
                 pressure='hPa',             $
                 temperature='K',            $
                 absunit='1/m',              $
                 plotfilename=plotfilename,  $
                 plotfileformat=2,           $
                 plotsum=2,                  $
                 plotyaxis='log',            $
                 avoid_tg=avoid_tg
;
goto, ende
;
; --------------------------------------------------------------------------

RTplot:
; T_B plots
plot_RT_1, jobname,                    $
           jobdir=jobdir,              $
           plotfilename=pabsname,      $
           plotyaxis='lin'
;
goto, ende
;
; --------------------------------------------------------------------------
;
TRAplot:
; transmission plots
plot_trans_1, jobname,               $
              jobdir=jobdir,         $
              plotfilename=pabsname, $
              yrange=[0.0,15.0]
;
;
goto, ende
;
; --------------------------------------------------------------------------
;
ende:
close, /all
END
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
