PRO sort_abs_2,absin,tgin,alt,absout,tgout,max_arr1,min_arr1
;+
;NAME: 
;      sort_abs_2.pro
;PURPOSE:
; rearrange the absorption and tag groups according to their
; absorption magnitude at altitude alt
;
; INPUT:
;     absin      : structure-absorption per tag group structure
;     tgin       : string array-names of tag groups
;     alt        : integer-index to altitude level
;
; OUTPUT:
;     absout     : structure-sorted absin
;     tgin       : string array-sorted tgin
;     max_arr1   : double array-maximum values found
;     min_arr1   : double array-minimum values found
;
; HISTORY:
;     2001-01-22 AvE created
;     2001-12-18 TKS adaption to plot_abs_per_tg_2
;-

;; dimension of structure
dim=(size(tgin,/DIMENSIONS))[0]

;; array to hold the max and min values of input data
max_arr=dblarr(dim)
min_arr=dblarr(dim)
;; array to hold the max and min values of output data
max_arr1=dblarr(dim)
min_arr1=dblarr(dim)

;; just copy the output arrays for the right dimensions
tgout=tgin
absout=absin

;; now fill the array with max values
for j=0,dim-1 do begin
    ;; we have to use this string exectution unfortunately, long story
    s1='max_arr['+string(j,format='(I0)')+'] = max(absin.mat'+$
      string(j,format='(I0)')+'[alt,*])'
    r1 = execute(s1)
    ;; min values
    s1='min_arr['+string(j,format='(I0)')+'] = min(absin.mat'+$
      string(j,format='(I0)')+'[alt,*])'
    r1 = execute(s1)
endfor

;; max_arr copy to overwrite
max1=max_arr

;; now sort the stuff 
for j=0,dim-1 do begin

    ;; find the index to max element
    a = max(max1,ind)

    ;; put the max value into the right position 
    tgout[j]=tgin[ind]
    s1='absout.mat'+string(j,format='(I0)')+' = absin.mat'+$
      string(ind,format='(I0)')
    r1 = execute(s1)

    ;; put the min and max values at right position
    min_arr1[j]=min_arr[ind]
    max_arr1[j]=max_arr[ind]

    ;; we have this value, set it to negative values
    max1[ind]=-1

endfor

end
