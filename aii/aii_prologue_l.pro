pro aii_prologue_l, color=color, ps=ps, prog, ls, thick, size=size, nomargin=nomargin
;----------------------------------------------------------------------
;   Standard Prologue for plotting
;----------------------------------------------------------------------
print,'Plotting '+prog+' ...'
print,'-----------------------------'
;print,'valid keyword parameters: /color /ps'
;print,'color=1: standard'
;print,'color=2: slide'
;print,'ps=1:    postscript'
;print,'ps=2:    encapsulated postscript'
;print,'defaults are: no color, x windows'
;print,''
;-----------------------------------------------------
;                 ls(i) is the linestyle	i > 0
;                 thick is the line thickness	
; 		  color 0 is black !
;-----------------------------------------------------
; History:
; SAB 23.01.98 Added keyword size = [xsize,ysize].
;              Default size is as it was.

thick      = 2.5   ; standard value of thick
slidethick = 5.    ; value of thick to use for slides

;; this messes up eps plots, thus allow to remove
if not keyword_set(nomargin) then begin
    !x.margin = [7.,2.]
    !y.margin = [3,2]
endif

if not keyword_set(size) then begin
    xsize=26
    ysize=19
endif else begin
    xsize=size[0]
    ysize=size[1]
endelse

if keyword_set(ps) then begin
    thick      = thick*2.
    slidethick = slidethick*2.
    set_plot,'ps',/interpolate
    !P.charsize = 2.
    if (ps NE 2) then begin 
        device,/landscape
        device,xoffset=1
        device,xsize=xsize
        device,ysize=ysize
        device,encapsulated=0
    end else begin
        device,/portrait
        device,xsize=xsize
        device,ysize=ysize
        device,xoffset=0
        device,/encapsulated
    end
    if keyword_set(color) then begin
        device,/color
        device,bits=8
        if (color EQ 1) then begin
            if (ps EQ 1) then 	device,filename=prog+'.cps'  $
            else 			device,filename=prog+'.ecps' 
        end else begin
            if (ps EQ 1) then 	device,filename=prog+'.fps'  $
            else			device,filename=prog+'.efps' 
        end
    end else begin
;        device,color=0
;; Make plot still colored!
        device,/color
        device,bits=8
        if (ps EQ 1) then		device,filename=prog+'.ps'   $
        else			device,filename=prog+'.eps'
    end
end else begin
    set_plot,'x'

    dummy         = !P.color            	; make background white
    !P.color      = !P.background
    !P.background = dummy

    !P.charsize = 1.12
    window,/free,xsize=fix(xsize*30),ysize=fix(ysize*30)
end

if keyword_set(color) then begin
    mycolor,slide = color - 1
    ls = intarr(31)
    ls = 0 * ls
    if (color EQ 2) then thick = slidethick
end else begin
    if keyword_set(ps) then begin
; Still use color!
        mycolor,slide = color - 1
    endif else begin
; Still use color!
        mycolor,slide = color - 1
    endelse

    lines = indgen(6)
    ls = [0,lines,lines,lines,lines,lines]
end

;; select font
if keyword_set(ps) then begin
    !p.font=0
    device,/helvetica
    device,bold = 0
    ; assign fonts !7 to symbol
    device, /SYMBOL, FONT_INDEX=7
    if keyword_set(color) then begin
        if color EQ 2 then begin
            device,/bold
        endif
    endif 
endif else begin
    !p.font = - 1
    if keyword_set(color) then begin
        if color EQ 2 then begin
            !P.FONT = 17
        endif
    endif 
endelse


;--------end-of-prologue---------------------------------------------
end
