PRO jpl_scaling,tag_arr,abun,index,count,for_species,hit_arr1,hit_arr2,hit_arr3,$
                scaled_abun
; checks whether the jpl abundancies must be scaled (maximum abundancy
; given in jpl catalogue eq 1). Performs the scaling with the maximum
; hitran abundancy found for that species. gives back the original jpl
; abundancy if the species in not present in hitran or no max
; abundancy with 1 was found in jpl cat

i = where(abun[index] eq 1,counter)

if counter gt 0 then begin

    j=where(hit_arr1 eq for_species,counter1)

    if counter1 gt 0 then fac=max(double(hit_arr3[j])) else fac = 1.0

    for k=0,count-1 do begin
        scaled_abun[k] = abun[index[k]] * fac
    endfor

endif else for k=0,count-1 do scaled_abun[k] = abun[index[k]]

end

        


PRO cat_abundancies,output=output

; reads a file which connects the jpl tag numbers with the species 
; name we use, this file is extracted from the jpl catalogue,
; currently placed at:
;
; /pool/lookup/jpl/cat7_00/abundancies/tag_species.jpl 
;
; opens then the corresponding jpl catalogue tag info
; file to read the abundancies, these files are extracted from the jpl
; web server
;
; keywords:
;        output    : string, if present writes output to the filename given,
;                    otherwise to screen


if n_elements(output) gt 0 then begin
    unit=3
    openw,3,output
    print,'Generating output file: ',output
endif else unit=-1

; array to hold tag numbers, jpl name, our name , resize array later
tag_arr=strarr(3,1000)


; everything with a # at the beginning is a comment, skip all of them
str='#'
openr,1,'/pool/lookup/jpl/cat7_00/abundancies/tag_species.jpl'
while strpos(str,'#') ge 0 do begin
    readf,1,str
end

; read the first entry into tag_arr
i=0
str=str_sep(strtrim(strcompress(str),2),' ')
tag_arr[0,i]=str[0] & tag_arr[1,i]=str[1] & tag_arr[2,i]=str[2]

; now read the complete file into tag_arr
while not eof(1) do begin
    str=''
    readf,1,str
    str=str_sep(strtrim(strcompress(str),2),' ')
    i=i+1
    tag_arr[0,i]=str[0] & tag_arr[1,i]=str[1] & tag_arr[2,i]=str[2]
end

close,1

; resize the tag_arr
tag_arr=tag_arr[*,0:i]

; generate the abundancy array
abun=dblarr(i)
for j=0,i-1 do abun[j]=100

; start reading the abundancies from d tag .cat files

for j=0,i-1 do begin
    str=''
    openr,1,'/pool/lookup/jpl/cat7_00/doc/d'+string(tag_arr[0,j],format='(I6.6)')+'.cat'
    while strpos(str,'Isotope Corr.') lt 0 do readf,1,str

    ; extract the isotope correction
    str=str_sep(str,'&')

    ; some of the isotope corrections are empty, check that first
    str=strtrim(strcompress(str[1]),2)

    ; calculate the abundancy and write into abun array
    if strlen(str) eq 0 then abun[j]=100 $
    else begin
        m=double(str)
        if m ge 0 then m=m*(-1)
        abun[j] = 10.0^m
    endelse

    close,1
endfor





; species available in the forward program, taken form glob_def.c
; file
; species 44: clono2 is specified in jpl as clno3, changed name to clono2


for_species=['H2O','CO2','O3','N2O','CO','CH4','O2','NO',$
             'SO2','NO2','NH3','HNO3','OH','HF','HCl','HBr',$    
             'HI','ClO','OCS','H2CO','HOCl','N2','HCN','CH3Cl',$  
             'H2O2','C2H2','C2H6','PH3','COF2','SF6','H2S','HCOOH',$
             'HO2','O','ClONO2','NO+','null','null','null','null',$
             'null','null','OClO','null','null','BrO','null','H2SO4',$ 
             'Cl2O2']


; initialize the hitran array, well well well, could have been
; structures 
hit_arr1=strarr(1000)
hit_arr2=intarr(1000)
hit_arr3=strarr(1000)
k=0

;; jpl catalogue files
jpl_cat=['/pool/lookup/jpl/newcat/cat00100.tag',$
         '/pool/lookup/jpl/newcat/cat00200.tag',$
         '/pool/lookup/jpl/newcat/cat00300.tag',$
         '/pool/lookup/jpl/newcat/cat00400.tag',$
         '/pool/lookup/jpl/newcat/cat00500.tag',$
         '/pool/lookup/jpl/newcat/cat00600.tag',$
         '/pool/lookup/jpl/newcat/cat00700.tag',$
         '/pool/lookup/jpl/newcat/cat00800.tag',$
         '/pool/lookup/jpl/newcat/cat00900.tag',$
         '/pool/lookup/jpl/newcat/cat01000.tag',$
         '/pool/lookup/jpl/newcat/cat01500.tag',$
         '/pool/lookup/jpl/newcat/cat02000.tag',$
         '/pool/lookup/jpl/newcat/cat03000.tag',$
         '/pool/lookup/jpl/newcat/cat05000.tag']
         

; cycle through all the species defined in the forward model, and
; check whether they are given in hitran and jpl catalogue
for j=0,n_elements(for_species)-1 do begin

    if for_species[j] ne 'null' then begin

        ; output the species:
        printf,unit,'Species: ',for_species[j]

        ; give the mytran tag
        printf,unit,'Mytran Tag: ',j+1


        ; get the degrees of freedom from jpl cataolgue

        ; find the array elements of the current forward species
        index=where(for_species[j] eq tag_arr[2,*],count)

        if count ne 0 then begin
            for l=0, n_elements(jpl_cat) - 1 do begin
                str=''
                openr,1,jpl_cat[l]
                repeat begin
                    readf,1,str
                endrep until (strtrim(strmid(str,45,6),2) eq tag_arr[0,index[0]]) or eof(1)
                close,1
                if (strtrim(strmid(str,45,6),2) eq tag_arr[0,index[0]]) then begin
                    ; extract the degrees of freedom
                    printf,unit,'Degrees of Freedom: ',strmid(str,29,2)
                    printf,unit,'Cat Line: ',str
                    printf,unit,''
                    goto, shit
                endif
            end
        endif

        shit:printf,unit,''

        ; read the hitran database, this file is just a copy of the hitran
        ; abundancies given in cdrom: hitran92/tables.3 file. 

        ; Comments at the beginning of the file are identified with a # sign
        str='#'
        openr,1,'hitran_abundancies.txt'
        while strpos(str,'#') ge 0 do begin
            readf,1,str
        end


        ; now find the species and hitran tag in the file
        str=str_sep(strtrim(strcompress(str),2),' ')
        while (str[0] ne for_species[j]) and $
          not eof(1) do begin
            str=''
            readf,1,str
            str=str_sep(strtrim(strcompress(str),2),' ')
        end

        ; species is not present in hitran
        if eof(1) then begin
            printf,unit,''
            printf,unit,'Missing in hitran'
            printf,unit,''
            close,1
        endif else begin

            ; species number in hitran
            i1=strpos(str[1],'(')
            i2=strpos(str[1],')')
            printf,unit,'Hitran Tag: ',fix(strmid(str[1],1,i2-1))
            printf,unit,''


            ; species name in hitran
            hit_species=str[0]

            ; species exists in hitran, now read all the abundancies, 
            ; and the hitran molecule identifier of this species
            str=''
            readf,1,str
            str=strtrim(strcompress(str),2)

            printf,unit,'Hitran species/tag/Abundancy for species: '
        
            ; assure that we are not yet in the next species abundancies
            while (strlen(str) gt 0) and (strpos(str,'(')) lt 0 do begin
                str=str_sep(strtrim(strcompress(str),2),' ')
                hit_arr1[k]=hit_species & hit_arr2[k]=str[0] & hit_arr3[k]=str[1]
                printf,unit, hit_arr2[k],hit_arr1[k],hit_arr3[k],format='(I8,A15,A15)'
                k=k+1
                if not eof(1) then begin
                    str=''
                    readf,1,str
                endif else str=''
            end
            printf,unit,''
            close,1
    
        endelse


        ; now do the jpl cat

        ; find the array elements of the current forward species
        index=where(for_species[j] eq tag_arr[2,*],count)

        ; does this species exists?
        if count eq 0 then begin
            printf,unit,'Missing in JPL catalogue'
            printf,unit,''
            printf,unit,'------------------------------------------------------------'
        endif else begin
            ; species exists in jpl cat, output all tag numbers, the original jpl
            ; name, the abundancy according to  jpl, and the scaled abundancy, if an
            ; abundancy of 1 is assumed in jpl cat. the scaling factor is the
            ; maximum abundancy of the hitran abundancies found for this species
            printf,unit,'JPL Tag numbers/name/Abundancy/scaled Abundancy for species'

            ; check whether scaling of jpl abundancies is necessary
            scaled_abun=dblarr(count)
            jpl_scaling,tag_arr,abun,index,count,for_species[j],hit_arr1,hit_arr2,hit_arr3,scaled_abun

            ; output the species, found and calculated abundancies
            for k=0,count-1 do begin
                printf,unit,tag_arr[0,index[k]],tag_arr[1,index[k]],$
                  abun[index[k]],scaled_abun[k],format='(A8,A15,G15.8,G15.8)'
            endfor
            printf,unit,''
            printf,unit,'------------------------------------------------------------'
            printf,unit,''
        endelse

    endif
    
end

; close the output file
if n_elements(output) gt 0 then close,3

end
