PRO jpl_scaling,tag_arr,index,count,for_species,hit_arr1,hit_arr2,hit_arr3
; checks whether the jpl isotopic ratios must be scaled (maximum isotopic ratio
; given in jpl catalogue eq 1). Performs the scaling with the maximum
; hitran isotopic ratio found for that species. gives back the original jpl
; isotopic ratio if the species in not present in hitran or no max
; isotopic ratio with 1 was found in jpl cat

i = where(tag_arr[index].tag_abun eq 1,counter)

if counter gt 0 then begin

    j=where(hit_arr1 eq for_species,counter1)

    if counter1 gt 0 then fac=max(double(hit_arr3[j])) else fac = 1.0

    for k=0,count-1 do begin
        tag_arr[index[k]].tag_abun_s = tag_arr[index[k]].tag_abun * fac
    endfor

endif else for k=0,count-1 do tag_arr[index[k]].tag_abun_s = tag_arr[index[k]].tag_abun

end

        


PRO cat_isotopic_ratio,output=output

; generates an output about found isotopic ratios from JPL and
; HITRAN. Several files are read in order to do this:
;
; tag_species.jpl : refers the jpl tag number to the ARTS names (only
;                   the names species are separated in the listing)
; hitran_isotopic_ratio.txt : isotopic ratios from hitran species
;
; jpl files c+tag number+.cat files for isotopic rations
; jpl files d+tag number+.cat files for degrees of freedom
;
; The directories where these files can be found are given below.
;
; For further info refer to the README file.
;
;
; keywords:
;        output    : string, if present writes output to the filename given,
;                    otherwise to screen


;; EDIT THIS:

;;---------------------------------------------------------------------------

;; where to find the jpl catalogue files
d_tag_files='/pool/lookup/jpl/cat7_00/doc/'  ; description jpl files
c_tag_files='/smiles_local/axel/'             ; line jpl files
         
; all species available in the forward program, taken from glob_def.c
; file + several other new ones.
; species 44: clono2 is specified in jpl as clno3, changed name to
; clono2
; add new species (in agreement with the tag_species.jpl file to this
; array. 


for_species=['H2O','CO2','O3','N2O','CO','CH4','O2','NO',$
             'SO2','NO2','NH3','HNO3','OH','HF','HCl','HBr',$    
             'HI','ClO','OCS','H2CO','HOCl','N2','HCN','CH3Cl',$  
             'H2O2','C2H2','C2H6','PH3','COF2','SF6','H2S','HCOOH',$
             'HO2','O','ClONO2','NO+','null','null','null','null',$
             'null','null','OClO','null','null','BrO','null','H2SO4',$ 
             'Cl2O2','HOBr','C2H4','OBrO','ClNO2','HOONO2','PO2','PS',$
            'C2S','C3O','HC3N','HNC3','HC2NC','C3N','NS','H2CS','PO','CS',$
            'PN','HNCO','HCP','CP','CH2CO','CH3CN','NH2CN','CH3C2H','SH',$
            'HNO','CH2NH','HCO','HNC','CN','C2H','NH','CH']

;;---------------------------------------------------------------------------


if n_elements(output) gt 0 then begin
    unit=3
    openw,3,output
    print,'Generating output file: ',output
endif else unit=-1

; structure to hold tag numbers, jpl name, our name, abundance,
; degrees of freedom
; initialize large, resize later
jpl_str={jpl_entries,$
         tag_nr      : ' ',$        ;; jpl tag number
         tag_name    : ' ',$        ;; jpl name to the tag number
         tag_arts    : ' ',$        ;; tag in arts name convention
         tag_abun    : 100.0,$      ;; abundance of the species
         tag_abun_s : 100.0,$      ;; scaled abundance of the species
         tag_dfr     : -1}          ;; degrees of freedom 

tag_arr = replicate({jpl_entries},1000)


;; read the jpl tag number, name, and the corresponding arts name

; everything with a # at the beginning is a comment, skip all of them
str='#'
openr,1,'tag_species.jpl'
while strpos(str,'#') ge 0 do begin
    readf,1,str
end

; read the first entry into tag_arr, currently present in str
i=0
str=str_sep(strtrim(strcompress(str),2),' ')
tag_arr[i].tag_nr=str[0] & tag_arr[i].tag_name=str[1] & tag_arr[i].tag_arts=str[2]

; now read the complete file into tag_arr
while not eof(1) do begin
    str=''
    readf,1,str
    str=str_sep(strtrim(strcompress(str),2),' ')
    i=i+1
    tag_arr[i].tag_nr=str[0] & tag_arr[i].tag_name=str[1] & tag_arr[i].tag_arts=str[2]
end
close,1

; resize the tag_arr
tag_arr=tag_arr[0:i]

; start reading the isotopic ratios from d tag .cat files
for j=0,i-1 do begin
    str=''
    openr,1,d_tag_files+'d'+string(tag_arr[j].tag_nr,format='(I6.6)')+'.cat'
    while strpos(str,'Isotope Corr.') lt 0 do readf,1,str

    ; extract the isotope correction
    str=str_sep(str,'&')

    ; some of the isotope corrections are empty, check that first
    str=strtrim(strcompress(str[1]),2)

    ; calculate the isotopic ratio and write into abun
    if strlen(str) eq 0 then tag_arr[j].tag_abun=100 $
    else begin
        m=double(str)
        if m ge 0 then m=m*(-1)
        tag_arr[j].tag_abun = 10.0^m
    endelse

    close,1
endfor


; start reading the degrees of freedom from c tag .cat files
; and perform a check with all data in this tag file, whether the dfr
; are identical
for j=0,i-1 do begin
    str=''
    openr,1,c_tag_files+'c'+string(tag_arr[j].tag_nr,format='(I6.6)')+'.cat'
    repeat begin
        readf,1,str
        tag_arr[j].tag_dfr=fix(strmid(str,29,2))
        chk=tag_arr[j].tag_dfr
        if not (tag_arr[j].tag_dfr eq chk) then begin
            print,'Error: Found degrees of freedom for one tag not identical.'
            print,'       Tag Nr: :'+tag_arr[j].tag_nr
            stop
        endif
    endrep until eof(1)
    close,1
endfor





; initialize the hitran array, well well well, could have been
; structures 
hit_arr1=strarr(1000)
hit_arr2=intarr(1000)
hit_arr3=strarr(1000)
k=0



; cycle through all the species defined in the forward model, and
; check whether they are given in hitran and jpl catalogue
for j=0,n_elements(for_species)-1 do begin

    if for_species[j] ne 'null' then begin

        ; output the species:
        printf,unit,'Species: ',for_species[j]

        ; print the degrees of freedom from jpl cataolgue
        ; find the array elements of the current forward species
        index=where(for_species[j] eq tag_arr[*].tag_arts,count)
        if count ne 0 then printf,unit,'Degrees of freedom: '+$
          string(tag_arr[index[0]].tag_dfr,format='(I2)')$
          else printf,unit,'Degrees of freedom: Not available in JPL'
        
        ; read the hitran database, this file is mainly a copy of the hitran
        ; isotopic ratios given in cdrom: hitran92/tables.3 file. 

        ; Comments at the beginning of the file are identified with a # sign
        str='#'
        openr,1,'hitran_isotopic_ratio.txt'
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

            ; species exists in hitran, now read all the isotopic ratios, 
            ; and the hitran molecule identifier of this species
            str=''
            readf,1,str
            str=strtrim(strcompress(str),2)

            printf,unit,'Hitran species/tag/Isotopic Ratio for species: '
        
            ; assure that we are not yet in the next species isotopic ratio
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
        index=where(for_species[j] eq tag_arr[*].tag_arts,count)

        ; does this species exists?
        if count eq 0 then begin
            printf,unit,'Missing in JPL catalogue'
            printf,unit,''
            printf,unit,'------------------------------------------------------------'
        endif else begin
            ; species exists in jpl cat, output all tag numbers, the
            ; original jpl name, the isotopic ratio according to jpl,
            ; and the scaled isotopic ratio, if an isotopic ratio of 1
            ; is assumed in jpl cat. the scaling factor is the maximum
            ; isotopic ratio of the hitran isotopic ratio found for this
            ; species
            printf,unit,'JPL Tag numbers/name/Isotopic Ratio/scaled Isotopic Ratio for species'

            ; check whether scaling of jpl isotopic ratios is necessary
            jpl_scaling,tag_arr,index,count,for_species[j],hit_arr1,hit_arr2,hit_arr3

            ; output the species, found and calculated isotopic ratios
            for k=0,count-1 do begin
                printf,unit,tag_arr[index[k]].tag_nr,tag_arr[index[k]].tag_name,$
                  tag_arr[index[k]].tag_abun,tag_arr[index[k]].tag_abun_s,format='(A8,A15,G15.8,G15.8)'
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
