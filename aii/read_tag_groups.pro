; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
function find_keyword , kw , cf , n
;;----------------------------------------------------------------------
;; Finds the keyword kw in array of strings cf of size n
;;----------------------------------------------------------------------
for i=0,n-1 do begin
   if (strpos(cf[i],kw) EQ 0) then return, i
endfor
return, -1
end


;; #########################################################################


PRO OPEN_READ_CONTROLFILE, controlfile, cf, nlines


print,'OPEN_READ_CONTROLFILE> Reading ARTS controlfile: ',controlfile


;; get the # of lines in the controlfile
spawn ,'wc -l '+controlfile, result
nlines = 1                              ; make nlines integer
reads, result, nlines
cf = strarr(nlines)

;; read the whole controlfile into array cf
openr,u,controlfile,/get_lun
readf,u,cf
close,u
free_lun,u

;; remove leading and trailing and inbetween blanks
cf = strtrim(strcompress(cf,/REMOVE_ALL),2)

;; remove everything behind a #
for i=0,nlines-1 do begin
    pos=strpos(cf[i],'#')
    if pos ge 0 then cf[i]=strmid(cf[i],0,pos) 
endfor

end


;; #########################################################################


PRO STRTOTALCOUNT, str, searchstr, cnt, pos

;; counts the number of occurences of searchstr in str

;; Start searching in character position 0.
I = 0

;; Number of occurrences found.
cnt = 0

;; save positions, there should not be more than 200.
pos=intarr(200)

;; Search for an occurrence.
WHILE (I NE -1) DO BEGIN 

    ;;If one is found, count it and advance to next character position.
    I = STRPOS(str, searchstr, I)

    IF (I NE -1) THEN BEGIN 

        ;; save position 
        pos[cnt]=I

        ;; Update counter.
        cnt = cnt + 1

        ;; Increment I so as not to count the same instance of searchstr twice.
        I = I + 1

    ENDIF
ENDWHILE

;; resize pos
if cnt ne 0 then pos=pos[0:cnt-1]

END


;; #########################################################################


PRO READ_TAG_GROUPS, controlfile, searchpattern, tg

;; opens an arts controlfile and returns the tag groups found after
;; the searchpattern in the string array tg. This version is not dummy
;; approved, but can handle comment lines, varies definition of tag
;; groups.
;;
;; Example: read_tag_groups,'example.arts','tag_groupsDefine',tg
;; 
;;          gives back all the tag groups defined in
;;          tag_groups_Define, the syntax expected is something like:
;;          tag_groupsDefine{ ["tag group1", "tag group2",
;;                             "tag group1" ] }
;;          where the [ ] brackets are used to search for all tag
;;          groups. 
;;
;; INPUT:
;;        controlfile    : String-Arts controlfile name (including .arts)
;;        searchpattern  : String-generally the name of a workspace method
;;                         that defines tag groups
;;
;; OUTPUT:
;;        tg             : string array-all found tag groups
;;
;;
;; HISTORY:
;; 2001-01-22 Ave Created.
;;

print,' '
print,'READ_TAG_GROUPS> controlfile  ='+controlfile
print,'READ_TAG_GROUPS> searchpattern='+searchpattern
print,' '

;; read controlfile into array cf, which holds nlines lines
open_read_controlfile, controlfile, cf, nlines

;; now search for the searchpattern
i = find_keyword(searchpattern,cf,nlines)
print,'READ_TAG_GROUPS> find_keyword: i=',i

;; nothing found?
if i lt 0 then begin
    print,'Search pattern: ',searchpattern,' not found in Controlfile: ',controlfile
    stop
endif

;; now search for the [] brackets
while strpos(cf[i],'[') lt 0 do i=i+1

;; ok, [ is found, now search for the closing one ] and write all
;; lines into a string array, initialize this array with 100 and
;; resize later
tg=strarr(100)
k=i

while strpos(cf[i],']') lt 0 do begin
    tg[i-k]=cf[i]
    i=i+1
endwhile 

;; put the last element into tg
tg[i-k] = cf[i]
;; total number of lines read
tnl=i-k

;; how many tag groups are out there? This is not equals tnl, since
;; there could be more than 1 tag group per line
tot_tg=0
for j=0,tnl do begin
    strtotalcount, tg[j], '"', cnt, pos
    tot_tg = tot_tg + cnt
endfor

;; total number of tag groups:
tot_tg = tot_tg/2


;; now allocate the right tg size and use temporary array tg1
tg1=strarr(tot_tg) 

;; remove everything outside of " "
tgi=0
for j=0,tnl do begin

    ;; get total number and position of " "
    strtotalcount, tg[j], '"', cnt, pos

    ;; we only have one tag group
    if cnt eq 2 then begin
        tg1[tgi] = strmid(tg[j],pos[0]+1,pos[1]-(pos[0]+1))
        tgi=tgi+1
    endif else begin
        if cnt gt 2 then begin
            ;; we have more than one tag group
            for k=0,cnt/2-1 do begin
                tg1[tgi]=strmid(tg[j],pos[k*2]+1,pos[k*2+1]-(pos[k*2]+1))
                tgi=tgi+1
            endfor
        endif
    endelse
endfor

;; copy the new tg groups back to tag
tg=tg1

end
;; #########################################################################
