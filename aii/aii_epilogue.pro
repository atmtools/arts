pro aii_epilogue, prog, ps, nostamp=nostamp, annotate=annotate
;-------------standard-epilogue-------------------------------
; The annotate option does not yet work correctly
;-------------------------------------------------------------

if (not keyword_set(nostamp)) then begin 
    if keyword_set(ps) then begin
;        if (ps EQ 2) then nostamp = 1 else nostamp = 0
;        Always plot the stamp
        nostamp = 0
    endif else nostamp = 0
endif

if (nostamp EQ 0) then begin
    ;; Get date:
    spawn,'date +"%y"',year
    spawn,'date +"%m"',month
    spawn,'date +"%d"',day
    datum = string(day)+'.'+string(month)+'.'+string(year)
    ;; Get user name:
    spawn,'whoami',who
    ;; Some defaults:
    case who[0] of
        'sbuehler' : who='SAB'
        'jo'       : who='JU'
        'axel'     : who='AvE'
        'bms'      : who='BMS'
        'cverdes'  : who='CV'
        'tkuhn'    : who='TKS'
        default    : who=who[0]
    endcase
    xyouts,0,0,/normal,size=.5,who+' '+datum+' / '+prog
endif

if keyword_set(annotate) then begin
    mydevice = !D.NAME
    set_plot,'x'                ;change to X but remember device
    ;;call annotate widget with default annotation file name
    annotate,load_file=prog+'.anot.dat' 
    set_plot,mydevice           ;set back output device
endif

if keyword_set(ps) then begin
    device,/close_file
    ;;spawn,'nice -19 lpr -h idl.'+strlowcase(devname)
    set_plot,'x'
endif

;; Make plot group readable
spawn,'chmod g+rw '+prog+'*ps 2>/dev/null'



;----------------------------------------------------------------
end
