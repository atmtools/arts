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
;       run          (string)  action selection, valid strings are
;                              'VMR'    VMR plot
;                              'ABS'    absorption plot
;                              'RT'     TB plot
;                              'TRANS'  transmission plot
; OUTPUTS:
;       Postscript file with absorption plot.
;
; EXTERNAL CALLS:
;       plot_vmr_per_tg
;       plot_abs_per_tg_2
;       plot_RT_1
;       plot_trans_1
;
; MODIFICATION HISTORY:
;       04-27-2001  TKS  alpha version created 
;       03-04-2003  TKS  update of plotting actions
;
; ==========================================================================

close, /all

; input control:
if not keyword_set(jobname) then begin
    print,'showme> !!! ERROR! no ARTS job name given!'
    print,'showme> !!! terminate without action.'
    goto, ende
endif

if not keyword_set(jobdir) then begin
    spawn,'pwd ', jobdir
endif

if not keyword_set(plotfilename) then begin
    plotfilename = jobname
endif

if not keyword_set(altitude) then begin
    selalt = 1.0
endif else begin
    selalt = altitude
endelse

if not keyword_set(avoid_tg) then begin
    avoid_tg = ''
endif

if not keyword_set(run) then begin
    print,'showme> !!! ERROR! no action specified!'
    print,'showme    possible keywords for run:'
    print,'showme>   1. VMR   <--->  VMR plot'
    print,'showme>   2. ABS   <--->  absorption per tag plot'
    print,'showme>   3. RT    <--->  T_B plot'
    print,'showme>   4. TRANS <--->  transmission plot'
    print,'showme> !!! terminate without action.'
    goto, ende
endif


print, ' '
print, '----------------------------------------------------------'
print, 'showme> action           : '+run
print, 'showme> arts jobname     : '+jobname
print, 'showme> arts jobdir      : '+jobdir
print, 'showme> output filename  : '+plotfilename
;
spawn, 'rm -f '+plotfilename+'.ps'
print, ' '
print, 'showme> selected altitude: ',selalt,' km'
print, '----------------------------------------------------------'
print, ' '

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (STRUPCASE(run) EQ 'VMR')   then goto, vmrplot
if (STRUPCASE(run) EQ 'ABS')   then goto, absplot
if (STRUPCASE(run) EQ 'RT')    then goto, RTplot
if (STRUPCASE(run) EQ 'TRANS') then goto, TRAplot

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

absplot:
pabsname = plotfilename+'_abs'
plot_abs_per_tg_2, jobname,                  $
                 altitude=selalt,            $
                 jobdir=jobdir,              $
                 pressure='hPa',             $
                 temperature='K',            $
                 absunit='1/m',              $
                 plotfilename=plotfilename,  $
                 plotfileformat=4,           $
                 plotsum=2,                  $
                 plotyaxis='log',            $
                 avoid_tg=avoid_tg
;
goto, ende

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RTplot:
; T_B plots
plot_RT_1, jobname,                    $
           jobdir=jobdir,              $
           plotfilename=pabsname,      $
           plotyaxis='lin'
;
goto, ende

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAplot:
; transmission plots
plot_trans_1, jobname,               $
              jobdir=jobdir,         $
              plotfilename=pabsname, $
              yrange=[0.0,15.0]
;
;
goto, ende

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ende:
close, /all
END
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
