; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;

Function aii_where_str,str1,str2,dim
;+
;NAME:
;      aii_where_str
;PURPOSE:
;      similar to where, but works for string array str1 and str2,
;      returns -1 if not found, otherwise the index of str1 where entries
;      of str2 were found
;USAGE: 
;      result = aii_where_str(str1,str2,dim)
;ARGUMENTS:
;      str1,str2: string arrays; occurrences of str2 elements in
;                 str1 are looked for
;      dim: variable that contains the number of occurrences 
;-

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

; ==========================================================================
