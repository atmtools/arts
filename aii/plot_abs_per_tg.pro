Function where_str,str1,str2,dim
;; similar to where, but works for string array str1 and str2
;;
;; returns -1 if not found, otherwise the index of str1 where entries
;; of str2 were found

;; size of str1, str2
ns1=(size(str1,/dime))[0]
ns2=(size(str2,/dime))[0]

;; allocate an int array, resize it later to dim
x=intarr(ns1)
dim=0

;; now loop the 2 arrays
for i=0,ns1-1 do begin
    for j=0,ns2-1 do begin
        if str1[i] eq str2[j] then begin
            x[dim]=i
            dim=dim+1
        endif
    endfor
endfor

;; now resize and return
if dim eq 0 then x=-1 else x=x[0:dim-1]
return, x
end


PRO plot_abs_per_tg, jobname, f, abs, alt, altitude=altitude,read=read,$
                     epsfile=epsfile,psfile=psfile,$
                     add_to_title=add_to_title,$
                     avoid_tg=avoid_tg,color=color,cm=cm,$
                     xrange=xrange,yrange=yrange

;; reads an arts abs_per_tg file, and plots the different tag group
;; absorption, the tag groups are identified by reading the arts
;; controlfile. procedure does not plot zero absorption of tag
;; groups.
;;
;; INPUT:
;;     jobname         : jobname of calculation
;; 
;; OUTPUT: (only needed to use the read keyword)
;;     f               : frequency grid
;;     abs             : altitude per tag group
;;     alt             : altitude grid
;;
;; KEYWORDS:
;;     altitude        : make abs plot near that altitude, 
;;                       default: 25 km
;;     read            : do not read the f, abs, alt again
;;     epsfile         : string-name of eps file to generate
;;     psfile          : string-name of ps file to generate
;;     add_to_title    : string-is added to default title
;;     avoid_tg        : array of strings-which tag groups to plot
;;     color           : make color plot (actually always color)
;;     cm              : use cm^-1 instead of GHz
;;     xrange          : array of double-plot range
;;     yrange          : array of double-plot range
;;
;; HISTORY:
;;     2001-01-25 : Created Ave
;;


;; read the tag groups from controlfile
read_tag_groups,jobname+'.arts','tag_groupsDefine',tg

;; total number of tag groups
tot_tg= (size(tg,/DIMENSIONS))[0]

;; now reads the absorption per tag, and the altitudes and
;; frequencies, only ascii arrays are allowed currently
if not keyword_set(read) then begin
    print,'Reading data'
    f=read_datafile(jobname+'.f_mono.aa', /optimize)
    abs=read_datafile(jobname+'.abs_per_tg.aa', /optimize)
    alt=read_datafile(jobname+'.z_abs.aa', /optimize)
endif

;; check that the tag groups and the absoprtion per tag group have the
;; same dimension
if n_tags(abs) ne tot_tg then begin
    print,'Error: absorption per tag group and tag goups have not'
    print,'       the same dimension'
    stop
endif

;; altitude if not given about 25 km
if not keyword_set(altitude) then begin
    altitude = (where(alt gt 25000))[0]
endif else begin
    altitude = (where(alt gt altitude))[0]
endelse

;; rearrange the absorption and tag groups according to their
;; absorption magnitude, and get max and min values for plot range
sort_abs,abs,tg,altitude,abs1,tg1,max_arr,min_arr

;; are certain tag groups selected or should all be plotted?
;; we use min_arr values of zero for avoiding
if keyword_set(avoid_tg) then begin
    avoid_tg_in = where_str(tg1,avoid_tg,dim)
    if dim eq 0 then begin
        print,'Error: No matching avoid_tg found.'
    endif else begin
        min_arr[avoid_tg_in] = 0
    endelse
endif

;; plot range
if keyword_set(cm) then begin
    unit=30.0
    unitstr='[cm!u-1!n]'
endif else begin
    unit=1.0
    unitstr='[GHz]'
endelse
if not keyword_set(yrange) then $
  yrange=[min( min_arr(where(min_arr ne 0)))/10.0,max(max_arr)*10.0]
if not keyword_set(xrange) then $
  xrange=[min(f/1E9/unit),max(f/1E9/unit)]

;; make 2 plots, one with the legend, the other one with the data

;; check for eps output
if keyword_set(epsfile) then begin
    if not arg_present(epsfile) then filename = jobname $ 
    else filename=epsfile
    ps=2
endif else begin
    if keyword_set(psfile) then begin
        if not arg_present(psfile) then filename = jobname $ 
          else filename=psfile
        ps=1
    endif else begin
        filename=jobname
        ps=0
    endelse
endelse 

;; check for color output? -> always use colors
if not keyword_set(color) then color=1

;; use prologue and epilogue to produce output
;; save settings
P_ini = !P
aii_prologue_l, color=color, ps=ps, filename, ls, thick

;; create color and lineshape arrays
colors=intarr(tot_tg)
ls=intarr(tot_tg)
for i=0,tot_tg-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    ls[i]     = i mod  5    ;; we have  5 different lineshapes 
endfor

;; data plot
!P.multi=[0,2,1]
!X.margin = [2,-10]
!Y.margin = [4,4]

;; add to title string
if not keyword_set(add_to_title) then add_to_title=''

;; plot in units of GHz or cm^-1
plot,f/1E9/unit,abs1.mat0[altitude,*],yrange=yrange,ytype=1,$
  title='Absorption: '$
  +string(alt[altitude]/1000.0,format='(F6.2)')+' km'+add_to_title,$
  xrange=xrange,xtitle='Frequency '+unitstr,$
  ytitle='Absorption [1/m]',$
  color=colors[0],linestyle=ls[0],thick=thick,$
  xstyle=1,ystyle=1,/nodata;;,charsize=0.7*!P.charsize

index=0 ;; count only the ones we plot
nozeros=where(min_arr ne 0)
for j=0,tot_tg-1 do begin
    ;; check that this absorption tag group is not zero
    dum=where(nozeros eq j,count)
    if count eq 1 then begin
        s1='oplot, f/1E9/unit, abs1.mat'+string(j,format='(I0)')+$
          '[altitude,*],linestyle='+string(ls[index],format='(I0)')+$
          ', color='+string(colors[index],format='(I0)')+', thick=2*thick'
        index=index+1
        r1 = execute(s1)
    endif else begin
        print,'Avoided or found absorption is zero, no plotting: ',tg1[j]
    endelse
endfor

;; legend plot
!X.margin = [-20,0]
!Y.margin = [4,4]
plot,[0,0],xstyle=4,ystyle=4,/nodata,/noerase,/noclip
klegend_d,tg1[where(min_arr ne 0)],/fill,box=0,$
  usersym=usersym,charsize=1.1*!P.charsize,spacing=3.2,$
  psym=intarr(index),colors=colors[0:index-1],$
  line=ls[0:index-1],position=[0.6,0.5+index*0.03], $
  thick=2*thick;,textcolors=colors

!P.multi=0

;; close file/plot
aii_epilogue, filename, ps, /nostamp
;; restore settings
P_ini = !P    

end
