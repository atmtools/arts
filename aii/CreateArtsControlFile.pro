; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
; ########################### INTERNAL FUNCTIONS ########################### 
;
; **************************************************************************
; Name:     addfields2vector
;
; Purpose:  adds one field to a 1D vector
;
; Inputs:   a             1D vector of m elements
;           n             number of additional fields
;
; Output:   a             1D vector of m+n elements
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION addfields2vector, a, n

dim = SIZE(a, /N_DIMENSIONS)

IF (dim EQ 1) THEN BEGIN
    b = N_ELEMENTS(a)+n
    b[0:N_ELEMENTS(a)-1] = a
    RETURN, b
ENDIF ELSE BEGIN
    print, 'addfields2vector> !!! WARNING array dimension is not one, dim=',dim
    print, 'addfields2vector> !!! >>> output array = input array'
ENDELSE

RETURN, a
END
;
; **************************************************************************
; Name:     list_arts_tag
;
; Purpose:  check if the arts tag exists. These tags need additional
;           information in the arts control file.
;
; Inputs:   tag           arts tag name
;
; Output:   ntagfound     number of arts tags found with the input
;                         arts tag name 
;
; History:  2002-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
PRO list_arts_tag, path=path, ilevel=ilevel
;
;
; ---- SET DEGREE OF COMMENTS ---------------------------------------
IF NOT KEYWORD_SET(ilevel) THEN ilevel = 0 
;
;
; ========== GET INFORMATION FROM species_data.cc FILE ==============
;
; ---- CHECKS FOR TEMPORARY FILE ------------------------------------
path  = STRCOMPRESS(path, /REMOVE_ALL) ; remove whitespace
if (STRMID(path, STRLEN(path)-1, 1) NE '/') THEN path=path+'/'
spawn, 'rm -f temp_arts_file.txt'
spawn, 'egrep -w NAME\|REC '+path+'src/species_data.cc >> temp_arts_file.txt'
;
;
; ---- OPEN FILE FOR READING ----------------------------------------
openr, unit, 'temp_arts_file.txt', error=err, /get_lun
IF err NE 0 THEN BEGIN
    print,'list_arts_tag> !!! ERROR: can not open file >>temp_arts_file.txt<<'
    print,'list_arts_tag> !!! >>> please check the directory and file names! STOP!'
    FREE_LUN, unit & stop
ENDIF
;
; ---- COMMON BLOCK FOR ARTS TAG INFO -------------------------------
COMMON ARTSTAGINFO, ARTSTAGINFO_NAMES, ARTSTAGINFO_ISO, $
                    ARTSTAGINFO_TGM, ARTSTAGINFO_MOD
;
; ---- READ DATA FILE LINE BY LINE ----------------------------------
datarow   = ' '
i = -1
j = -1
NAMES = strarr(100)    ; size is a first guess
ISO   = strarr(100,25) ; size is a first guess
WHILE NOT eof(unit) DO BEGIN
    readf, unit, datarow
    p = str_sep(datarow, '"')   ; divides datarow into pieces
;   find lines with "NAME" in it
    pos1 = STRPOS(p[0], 'NAME')
;   find lines with "REC" in it
    pos2 = STRPOS(p[0], 'REC')
    IF ((pos1 GT 0) AND (N_ELEMENTS(p) GE 2) AND (pos2 LT 0)) THEN BEGIN ; new species
        i = i + 1 ; increase the species index
        j = -1 ; set the isotope index back
        NAMES[i] = p[1]
    ENDIF
    IF ((pos1 LT 0) AND (pos2 GT 0) AND (N_ELEMENTS(p) GE 2)) THEN BEGIN ; new isotope
        j = j + 1
        ISO[i,j] = p[1]
    ENDIF
ENDWHILE
FREE_LUN, unit
;
; ---- RESIZE THE ARRAYS --------------------------------------------
NAMES = NAMES[0:i]
ISO   = ISO[0:i,0:N_ELEMENTS(ISO[0,*])-1]
;
; ---- FILL COMMON BLOCK VARIABLES ----------------------------------
ARTSTAGINFO_NAMES = strarr(N_ELEMENTS(NAMES[*]))    ; size is a first guess
ARTSTAGINFO_ISO   = strarr(N_ELEMENTS(NAMES[*]), N_ELEMENTS(ISO[0,*])) ; size is a first guess
i = -1
j = -1
jmax = 0
FOR a = 0,N_ELEMENTS(NAMES[*])-1 DO BEGIN
    IF (STRLEN(ISO[a,0]) GT 0) THEN BEGIN 
        i = i + 1
        j = -1
        ARTSTAGINFO_NAMES[i] = NAMES[a]
        FOR b = 0,N_ELEMENTS(ISO[0,*])-1 DO BEGIN
            IF (STRLEN(ISO[a,b]) GT 0) THEN BEGIN 
                j = j + 1
                ARTSTAGINFO_ISO[i,j] = ISO[a,b]
                IF (j GT jmax) THEN jmax = j
            ENDIF
        ENDFOR
    ENDIF
ENDFOR
;
; ---- PRINT FOR CONTROL --------------------------------------------
;goto, list_arts_tag_jump_print 
ARTSTAGINFO_NAMES = ARTSTAGINFO_NAMES[0:i]
ARTSTAGINFO_ISO   = ARTSTAGINFO_ISO[0:i,0:jmax]
FOR i = 0,N_ELEMENTS(ARTSTAGINFO_NAMES)-1 DO BEGIN
    IF (ilevel GT 0) THEN print, 'list_arts_tag> ',i,': arts spec name=',ARTSTAGINFO_NAMES[i]
    FOR j = 0,N_ELEMENTS(ARTSTAGINFO_ISO[i,*])-1 DO BEGIN
        IF (STRLEN(ARTSTAGINFO_ISO[i,j]) GT 0) THEN $
           IF (ilevel GT 0) THEN print, 'list_arts_tag>    ',j,': iso=',ARTSTAGINFO_ISO[i,j]
    ENDFOR
ENDFOR
list_arts_tag_jump_print:
;
; ============ GET INFORMATION FROM CONTINUA.CC FILE ================
;
; ---- CHECKS FOR TEMPORARY FILE ------------------------------------
path  = STRCOMPRESS(path, /REMOVE_ALL) ; remove whitespace
if (STRMID(path, STRLEN(path)-1, 1) NE '/') THEN path=path+'/'
spawn, 'rm -f temp_arts_file.txt'
spawn, 'egrep CONTAGMODINFO '+path+'src/continua.cc >> temp_arts_file.txt'
;
; ---- OPEN FILE FOR READING ----------------------------------------
openr, unit, 'temp_arts_file.txt', error=err, /get_lun
IF err NE 0 THEN BEGIN
    print,'list_arts_tag> !!! ERROR: (2) can not open file >>temp_arts_file.txt<<'
    print,'list_arts_tag> !!! >>> please check the directory and file names! STOP!'
    FREE_LUN, unit & stop
ENDIF
;
; ---- READ DATA FILE LINE BY LINE ----------------------------------
datarow   = ' '
i = -1
jmax = 0
TAGS   = strarr(100)    ; size is a first guess
MODELS = strarr(100,10) ; size is a first guess
WHILE NOT eof(unit) DO BEGIN
    readf, unit, datarow
    datarow = STRCOMPRESS(datarow, /REMOVE_ALL)
    p = str_sep(datarow, ':')   ; divides datarow into pieces
    IF (N_ELEMENTS(p) EQ 2) THEN BEGIN
;       find tag name
        i = i + 1
        TAGS[i] = STRMID(p[0], 1+RSTRPOS(p[0], '*'))
;       find associated valid model variables
        p2 = str_sep(p[1], ',') ; divides datarow into pieces
        IF (N_ELEMENTS(p2) GT 0) THEN BEGIN
            FOR j = 0, N_ELEMENTS(p2)-1 DO MODELS[i,j] = p2[j]
        ENDIF
        IF (N_ELEMENTS(p2) GT jmax) THEN jmax = N_ELEMENTS(p2)
    ENDIF
ENDWHILE
FREE_LUN, unit
;
; ---- RESIZE THE ARRAYS --------------------------------------------
TAGS   = TAGS[0:i]
MODELS = MODELS[0:N_ELEMENTS(TAGS)-1,0:jmax]
;
; ---- PRINT FOR CONTROL --------------------------------------------
;goto, list_arts_tag_jump2_print 
FOR i = 0,N_ELEMENTS(TAGS)-1 DO BEGIN
    IF (ilevel GT 0) THEN print, 'list_arts_tag> ',i,': mod tag name=',TAGS[i]
    FOR j = 0,N_ELEMENTS(MODELS[i,*])-1 DO BEGIN
        IF (STRLEN(MODELS[i,j]) GT 0) THEN $
           IF (ilevel GT 0) THEN print, 'list_arts_tag>    ',j,': mod=',MODELS[i,j]
    ENDFOR
ENDFOR
list_arts_tag_jump2_print:
;
; ---- FILL COMMON BLOCK VARIABLES ----------------------------------
ARTSTAGINFO_TGM = TAGS
ARTSTAGINFO_MOD = MODELS
;
; ---- END OF PROCEDURE ---------------------------------------------
list_arts_tag_end:
END
;
; **************************************************************************
; Name:     ckeck_arts_tag
;
; Purpose:  check if the arts tag exists
;
; Inputs:   tag           arts tag name
;           model         tag associated model variable
;
; Output:   ntagfound     number of arts tags found with the input
;                         arts tag name 
;
; History:  2002-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION ckeck_arts_tag, tag, tagmodel
;
;print, 'ckeck_arts_tag> tag: >>',tag,'<<, model: >>',tagmodel,'<<'
;
; ---- ALL ARTS DEFINED TAGS -----------------------
; common block "ARTSTAGINFO" variables are: 
;   (a) ARTSTAGINFO_NAMES = arts supported species names
;   (b) ARTSTAGINFO_ISO   = arts supported species isotopes
;   (c) ARTSTAGINFO_TGM   = arts tag names wich have model parameter
;                           information to be specified within the 
;                           arts method cont_descriptionAppend
;   (d) ARTSTAGINFO_MOD   = valid parameter values of variable model
;                           of the arts method cont_descriptionAppend 
; This common block is fille by procedure 'list_arts_tag'
COMMON ARTSTAGINFO, ARTSTAGINFO_NAMES, ARTSTAGINFO_ISO, $
                    ARTSTAGINFO_TGM, ARTSTAGINFO_MOD
;
; --- TAG NAME SETTING -----------------------------
IF (STRLEN(tagmodel) GT 0) THEN tagmodel = STRCOMPRESS(tagmodel, /REMOVE_ALL)
tag       = STRCOMPRESS(tag, /REMOVE_ALL)
tagname   = str_sep(tag, '-')
;
; --- TAG NAME=SPECIES NAME ------------------------
IF (N_ELEMENTS(tagname) EQ 1) THEN BEGIN
    ntagfound = -1 ; at the end we have to add 1 to get the # of found tags
    taglist = intarr(1)
    FOR i = 0, N_ELEMENTS(ARTSTAGINFO_NAMES)-1 DO BEGIN
        IF (tagname[0] EQ ARTSTAGINFO_NAMES[i]) THEN BEGIN 
            ntagfound = ntagfound + 1
;            IF (ntagfound-1 GT N_ELEMENTS(taglist)) THEN taglist=addfields2vector(taglist,1)
;            taglist[ntagfound] = i
        ENDIF
    ENDFOR
;
;   ---- MODEL VERIFICATION ------------------------
    IF (( STRLEN(tagmodel) GT 0) AND (ntagfound GT 0)) THEN BEGIN
        print, 'ckeck_arts_tag> no model verification possible for tag >>',tag,'<<'
        print, 'ckeck_arts_tag> since not isotope name given.'
    ENDIF
    RETURN, ntagfound+1 ; number of tags found with this name
ENDIF 
;
; --- TAG NAME=SPECIES NAME AND ISOTOPE NAME -------
IF (N_ELEMENTS(tagname) EQ 2) THEN BEGIN
    ntagfound = 0
    n         = -1
    taglist_a = intarr(1)
    taglist_b = intarr(1)
    FOR i = 0, N_ELEMENTS(ARTSTAGINFO_NAMES)-1 DO BEGIN
        IF (tagname[0] EQ ARTSTAGINFO_NAMES[i]) THEN BEGIN 
            FOR j = 0, N_ELEMENTS(ARTSTAGINFO_ISO[i,*])-1 DO BEGIN
                IF (tagname[1] EQ ARTSTAGINFO_ISO[i,j]) THEN BEGIN 
;                    print, 'spec=',tagname[0],', iso=',tagname[1],$
;                           ' (',ARTSTAGINFO_NAMES[i],'-',ARTSTAGINFO_ISO[i,j],')'
                    n = n + 1
                    IF (n GT N_ELEMENTS(taglist_a)-1) THEN $
                       taglist_a=addfields2vector(taglist_a,(n-N_ELEMENTS(taglist_a)))
                    IF (n GT N_ELEMENTS(taglist_b)-1) THEN $
                       taglist_b=addfields2vector(taglist_b,(n-N_ELEMENTS(taglist_a)))
                    taglist_a[n] = i
                    taglist_b[n] = j
                ENDIF
            ENDFOR
        ENDIF
    ENDFOR
    ntagfound = n+1 ; number of found tags
;   ---- MODEL VERIFICATION ------------------------
    IF ( (STRLEN(tagmodel) GT 0) AND (ntagfound GT 0)) THEN BEGIN
        ntagmodfound = 0
        FOR n = 0, ntagfound-1 DO BEGIN ; loop over all resolved tags
            targettag = $
            ARTSTAGINFO_NAMES[taglist_a[n]]+'-'+ARTSTAGINFO_ISO[taglist_a[n],taglist_b[n]]
            IF (targettag NE tag) THEN BEGIN
                print, 'ckeck_arts_tag> !!! ERROR: wrong tag name found!'
                print, 'ckeck_arts_tag> input tag=',tag,', found tag=',targettag,', STOP!'
                STOP
            ENDIF 
;           loop over all tags with special model parameter values
            FOR i = 0, N_ELEMENTS(ARTSTAGINFO_TGM)-1 DO BEGIN
                IF (tag EQ ARTSTAGINFO_TGM[i]) THEN BEGIN
;                   loop over all species specific model parameter values
                    FOR j = 0, N_ELEMENTS(ARTSTAGINFO_MOD[i,*])-1 DO BEGIN
                        IF (tagmodel EQ ARTSTAGINFO_MOD[i,j]) THEN BEGIN 
                            ntagmodfound = ntagmodfound + 1
                        ENDIF
                    ENDFOR
                ENDIF
            ENDFOR
        ENDFOR
        RETURN, ntagmodfound  ; number of tags found with this name
    ENDIF 
    RETURN, ntagfound         ; number of tags found with this name
ENDIF
;
; ---- END -----------------------------------------
RETURN, 0 ; default number of tags found with this name
END
;
; **************************************************************************
; Name:     write_arts_control_file_tag
;
; Purpose:  writes appropriate arts tag information into a arts control file
;
; Inputs:   tags          vector of arts tag names
;
; Output:   flag          0=ok, 1=error occured
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION write_arts_control_file_tag, unit, tags
;
; ---- SET OUTPUT FLAG -----------------------------
flag = 1
;
; ---- CHECK ALL TAGS ------------------------------
ntags = N_ELEMENTS(tags)
FOR i = 0,ntags-1 DO BEGIN
    IF (ckeck_arts_tag(tags[i], '') NE 1) THEN BEGIN
        print, 'write_arts_control_file_tag> !!! ERROR, can not find arts tag=',tags[i],'!'
        print, 'write_arts_control_file_tag> !!! terminate now'
        STOP
    ENDIF
ENDFOR
;
; ---- CHECK IO ------------------------------------
IF (check_output_unit(unit) NE 0) Then BEGIN
    print, 'write_arts_control_file_tag> !!! ERROR, can not write to output file'
    print, 'write_arts_control_file_tag> !!! terminate now'
    STOP
ENDIF
;
; ---- WRITE TAGS INTO OUTPUT FILE -----------------
printf, unit,'tgsDefine{['
FOR i = 0,N_ELEMENTS(tags)-2 DO printf, unit, '"'+tags[i]+'", '
printf, unit, '"'+tags[N_ELEMENTS(tags)-1]+'"'
printf, unit,']}'
printf, unit, FORMAT='(A1)','#'
flag = 0
;
; ---- END OF FUNCTION -----------------------------
write_arts_control_file_tag_end:
RETURN, flag
END
;
;
; **************************************************************************
; Name:     write_arts_control_file_contdescription
;
; Purpose:  writes appropriate arts tag information into a arts control file
;
; Inputs:   tags          vector of arts tag names
;
; Output:   flag          0=ok, 1=error occured
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION write_arts_control_file_contdescription, unit, tags, $
                                                  models, $
                                                  userparameters
;
; ---- SET OUTPUT FLAG -----------------------------
flag = 1
;
; ---- CHECK IO ------------------------------------
IF (check_output_unit(unit) NE 0) Then BEGIN
    print, 'write_arts_control_file_contdescription> !!! ERROR, can not write to output file'
    print, 'write_arts_control_file_contdescription> !!! terminate now'
    STOP
ENDIF
;
; ---- ALL ARTS DEFINED TAGS -----------------------
; common block "ARTSTAGINFO" variables are: 
;   (a) ARTSTAGINFO_NAMES = arts supported species names
;   (b) ARTSTAGINFO_ISO   = arts supported species isotopes
;   (c) ARTSTAGINFO_TGM   = arts tag names wich have model parameter
;                           information to be specified within the 
;                           arts method cont_descriptionAppend
;   (d) ARTSTAGINFO_MOD   = valid parameter values of variable model
;                           of the arts method cont_descriptionAppend 
; This common block is fille by procedure 'list_arts_tag'
COMMON ARTSTAGINFO, ARTSTAGINFO_NAMES, ARTSTAGINFO_ISO, $
                    ARTSTAGINFO_TGM, ARTSTAGINFO_MOD
;
; ---- CHECK ALL TAGS ------------------------------
single = 0
ntags = N_ELEMENTS(tags)
FOR i = 0,ntags-1 DO BEGIN
;   make sure that this is a tag with an explicite isotope name
;   becaue only such specific taggs need additional information
;   via the arts method cont_descriptionAppend
    IF (STRPOS(tags[i], '-') GT 0) THEN BEGIN 
        flag = ckeck_arts_tag(tags[i], models[i])
;
;   ---- IF ONE TRUE ARTS TAG FOUND ----------------------
        IF (flag EQ 1) THEN BEGIN
;
;       ---- CHECK THAT A MODEL IS GIVEN IF NECESSARY ----
            foundmodel = 0
            FOR j = 0,N_ELEMENTS(ARTSTAGINFO_TGM)-1 DO BEGIN
                IF (ARTSTAGINFO_TGM[j] EQ tags[i]) THEN BEGIN
                    foundmodel = 1
                    print, 'write_arts_control_file_contdescription> ', $
                           'arts tag=',ARTSTAGINFO_TGM[j],', input tag=',tags[i]
                    IF (STRLEN(models[i]) LE 0) THEN BEGIN
                        print, 'write_arts_control_file_contdescription> !!! ERROR, in model detected!'
                        print, 'write_arts_control_file_contdescription> !!! this tag needs a model!'
                        print, 'write_arts_control_file_contdescription> !!! tag:   ',tags[i]
                        print, 'write_arts_control_file_contdescription> !!! model: ',models[i]
                        print, 'write_arts_control_file_contdescription> valid models are:',$
                          ARTSTAGINFO_MOD[j,*]
                        STOP
                    ENDIF
                ENDIF
            ENDFOR
;
;       ---- WRITE TAG INFO INTO OUTPUT FILE -------------
            if (single EQ 0) THEN BEGIN
                printf, unit,'cont_descriptionInit{}'
                single = 1      ; just write it once
            ENDIF
            printf, unit, 'cont_descriptionAppend{'
            printf, unit, '    tagname        = "'+tags[i]+'"'
            printf, unit, '    model          = "'+models[i]+'"'
            IF (models[i] EQ 'user') THEN BEGIN
                up = userparameters[*,i]
                ps = string(up) ; array of strings
                p = ' '
                FOR n = 0, N_ELEMENTS(ps)-1 DO BEGIN
                    s = STRCOMPRESS(ps[n], /REMOVE_ALL)
                    IF (STRLEN(s) GT 0) THEN p = p+','+s
                ENDFOR
                p = STRMID(p, 2)
                printf, unit, '    userparameters = [ '+p+' ]'
            ENDIF ELSE BEGIN
                printf, unit, '    userparameters = [ ]'
            ENDELSE
            printf, unit,'}'
        ENDIF ELSE BEGIN
            print, 'write_arts_control_file_contdescription> !!! ERROR in tag name check detected!'
            print, 'write_arts_control_file_contdescription> !!! tag:   ',tags[i]
            print, 'write_arts_control_file_contdescription> !!! model: ',models[i]
            STOP
        ENDELSE
    ENDIF
ENDFOR
printf, unit, FORMAT='(A1)','#'
flag = 0
;
; ---- END OF FUNCTION -----------------------------
write_arts_control_file_contdescription_end:
RETURN, flag
END
;
; **************************************************************************
; Name:     write_arts_control_file_lineshape
;
; Purpose:  writes appropriate arts tag information into a arts control file
;
; Inputs:   unit          output file unit
;           tags          vector of arts tag names
;           lineshapes    array of line shape information per tag
;
; Output:   flag          0=ok, 1=error occured
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION write_arts_control_file_lineshape, unit, tags, lineshapes
;
; ---- SET OUTPUT FLAG -----------------------------
flag = 1
;
; ---- VALID ARGUMENTS FOR LINE SHAPE NAMES --------
; these vectors have to be updated if any change is made
; in the arts source code since theses vectors are not automatically 
; generated from the actual source code.
valid_line_shapes_n = ['no_shape', $
                       'Doppler',  $
                       'Lorentz',  $
                       'Voigt_Kuntz3', 'Voigt_Kuntz4', 'Voigt_Kuntz6', $
                       'Voigt_Drayson', $
                       'Rosenkranz_Voigt_Drayson', 'Rosenkranz_Voigt_Kuntz6']
valid_line_shapes_f = ['no_norm', 'linear', 'quadratic']
;
; ---- CHECKS --------------------------------------
IF (N_ELEMENTS(lineshapes[0,*]) NE N_ELEMENTS(tags)) THEN BEGIN
    print, ' write_arts_control_file_lineshape> !!! ERROR in line shape specification!'
    print, ' write_arts_control_file_lineshape> !!! inconsistency in vector size of input tag vector '
    print, ' write_arts_control_file_lineshape> !!! and input line shape array '
    print, ' write_arts_control_file_lineshape> !!! line shape array size: ',$
           size(lineshapes, /N_DIMENSIONS)
    print, ' write_arts_control_file_lineshape> !!! tag vector size: ',$
           size(tags, /N_DIMENSIONS)
    RETURN, flag
ENDIF
;
FOR i = 0, N_ELEMENTS(lineshapes[0,*])-1 DO BEGIN
    ok = 1
    FOR j = 0, N_ELEMENTS(valid_line_shapes_n)-1 DO $
      IF (lineshapes[0,i] EQ  valid_line_shapes_n[j]) THEN ok = 0
    IF (ok NE 0) THEN BEGIN
        print, ' write_arts_control_file_lineshape> !!! ERROR in line shape specification!'
        print, ' write_arts_control_file_lineshape> !!! input tag:',$
        tags[i],', input line shape:',lineshapes[0,i]
        RETURN, flag
    ENDIF
ENDFOR
;
; ---- WRITE INTO CONTROL FILE ---------------------
printf, unit, 'lineshape_per_tgDefine{'
printf, unit, '    shape               = ['
FOR i = 0, N_ELEMENTS(lineshapes[0,*])-1 DO BEGIN
    IF (i LT N_ELEMENTS(lineshapes[0,*])-1) THEN a=',' ELSE a=''
    printf, unit, '                           "'+lineshapes[0,i]+'"'+a
ENDFOR
printf, unit, '                          ]'
printf, unit, '    normalizationfactor = ['
FOR i = 0, N_ELEMENTS(lineshapes[1,*])-1 DO BEGIN
    IF (i LT N_ELEMENTS(lineshapes[1,*])-1) THEN a=',' ELSE a=''
    printf, unit, '                           "'+lineshapes[1,i]+'"'+a
ENDFOR
printf, unit, '                          ]'
printf, unit, '    cutoff              = ['
FOR i = 0, N_ELEMENTS(lineshapes[2,*])-1 DO BEGIN
    IF (i LT N_ELEMENTS(lineshapes[2,*])-1) THEN a=',' ELSE a=''
;    reads, lineshapes[2,i], b
    printf, unit, '                           '+lineshapes[2,i]+a
ENDFOR
printf, unit, '                          ]'
printf, unit, '}'
;
;
; ---- SETOUTPUT FLAG ------------------------------
flag = 0
;
; ---- END OF FUNCTION -----------------------------
write_arts_control_file_lineshape_end:
RETURN, flag
END
;
; **************************************************************************
; Name:     write_arts_control_file_vmr
;
; Purpose:  writes appropriate VMR information into the arts control
;           file with the arts method raw_vmrsReadFromFiles. 
;
; Inputs:   unit            output file unit
;           tags            vector of arts tag names
;           vmrtagnames     array of tag names
;           vmrfilenames    array of VMR file names associated with vmrtagnames
;           vmrbasename     base name for VMR files. important for the
;                           tags which are in the variable tags but not in the
;                           variable vmrtagnames specified
;
; Output:   flag            0=ok, 1=error occured
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION write_arts_control_file_vmr, unit, tags, $
                                      vmrtagnames, vmrfilenames, vmrbasename 
;
; ---- SET OUTPUT FLAG -----------------------------
flag = 1
;
; ---- CHECKS --------------------------------------
IF (N_ELEMENTS(vmrtagnames) NE N_ELEMENTS(vmrfilenames)) THEN BEGIN
    print,'write_arts_control_file_vmr> ERROR: the two vectors vmrtagnames and'
    print,'write_arts_control_file_vmr> vmrfilenames have not the same size!'
    print,'write_arts_control_file_vmr> size of vmrtagnames:',SIZE(vmrtagnames, /N_DIMENSIONS)
    print,'write_arts_control_file_vmr> size of vmrfilenames:',SIZE(vmrfilenames, /N_DIMENSIONS)
    print,'write_arts_control_file_vmr> STOP here!'
    STOP
ENDIF
IF (N_ELEMENTS(vmrtagnames) GT N_ELEMENTS(tags)) THEN BEGIN
    print,'write_arts_control_file_vmr> ERROR: the vector vmrtagnames!'
    print,'write_arts_control_file_vmr> Size of vmrtagnames is > than the size of vector tags!'
    print,'write_arts_control_file_vmr> size of tags       :',SIZE(tags, /N_DIMENSIONS)
    print,'write_arts_control_file_vmr> size of vmrtagnames:',SIZE(vmrtagnames, /N_DIMENSIONS)
    print,'write_arts_control_file_vmr> STOP here!'
    STOP
ENDIF
n = 0
FOR i = 0, N_ELEMENTS(vmrtagnames)-1 DO BEGIN
    FOR j = 0, N_ELEMENTS(tags)-1 DO $
      IF (vmrtagnames[i] EQ tags[j]) THEN n = n+1  
ENDFOR
IF (n NE N_ELEMENTS(vmrtagnames)) THEN BEGIN
    print,'write_arts_control_file_vmr> ERROR: inconsistency in the vector vmrtagnames!'
    print,'write_arts_control_file_vmr> Not all the tags stated in vmrtagnames are also'
    print,'write_arts_control_file_vmr> found in the vector tags. STOP here!'
    STOP
ENDIF
FOR i = 0, N_ELEMENTS(vmrfilenames)-1 DO BEGIN
    IF (N_ELEMENTS(check_input_file(vmrfilenames[i],'')) NE 1) THEN BEGIN
        print,'write_arts_control_file_vmr> ERROR: unresolved VMR file names in variable vmrfilenames!'
        print,'write_arts_control_file_vmr> For ',vmrfilenames[i],' there are ',$
              check_input_file(vmrfilenames[i],''),' file(s) found. STOP here!'
        STOP
    ENDIF
ENDFOR
;
; ---- WRITE INTO CONTROL FILE ---------------------
printf, unit, 'raw_vmrsReadFromFiles'
printf, unit, '  {seltags   = ['
FOR i = 0, N_ELEMENTS(vmrtagnames)-1 DO BEGIN 
    IF (i LT N_ELEMENTS(vmrtagnames)-1) THEN a=',' ELSE a=''
    printf, unit, '                 "'+vmrtagnames[i]+'"'+a
ENDFOR
printf, unit, '               ]'
printf, unit, '   filenames = ['
FOR i = 0, N_ELEMENTS(vmrfilenames)-1 DO BEGIN 
    IF (i LT N_ELEMENTS(vmrfilenames)-1) THEN a=',' ELSE a=''
    printf, unit, '                 "'+vmrfilenames[i]+'"'+a
ENDFOR
printf, unit, '               ]'
printf, unit, '   basename  =  "'+vmrbasename+'"'
printf, unit, '  }'
;
; ---- SET OUTPUT FLAG -----------------------------
flag = 0
;
; ---- END OF FUNCTION -----------------------------
write_arts_control_file_vmr_end:
RETURN, flag
END
;
;
;############################# MAIN FUNCTION ################################
;
;
; **************************************************************************
; Name:     CreateArtsControlFile
;
; Purpose:  the bunch of procedures and functions built up an arts
;           control file and makes simultaneously some consistency checks.
;
; Inputs:   buffgastag    buffer gas arts tag name vector
; Inputs:   datafile            full path/file name
;           tags                full arts tag name and add. information
;                               e.g. 'tagname.model'
;           tag_models          tag associated model
;           tag_userparameters  tag associated userparameters (only for model='user')
;           frange              frequency range for fit [Hz]
;           debug               if selected turns additional printed
;                               information on (directed to stddev)
;           artsjobpath         subdirectory where the user runs the
;                               arts job
;           artspath            arts subdirectory path
;           catname             spectroscopic line catalog file name
;           catformat           spectroscopic line catalog format
;           catfmin             spectroscopic line catalog lower frequency limit [Hz]
;           catfmax             spectroscopic line catalog upper frequency limit [Hz]
;           lineshapes          line shape specifications (shape,
;                               normalization factor, and cutoff)
;                               lineshapes = strarr(3, N_ELEMENTS(tags))
;           frangemin           lower frequency limit of the frequency grid [Hz]
;           frangemax           upper frequency limit of the frequency grid [Hz]
;           frangesteps         number of frequency grid points
;           frangefile          input frequency gride file 
;           ptzfile             input file with
;                               pressure-temperature-altitude information
;                               in umits of [Pa], [K], and [m]
;           vmrtagnames         tag names used for variable seltags in arts method
;                               raw_vmrsReadFromFiles 
;           vmrfilenames        VMR file names used for variable filenames in arts method
;                               raw_vmrsReadFromFiles 
;           vmrbasename         base name used for variable basename in arts method
;                               raw_vmrsReadFromFiles 
;           vmrscenarioname     base name used for variable basename in arts method
;                               raw_vmrsReadFromScenario
;           prangemin           lower pressure limit of the pressure grid [Pa]
;           prangemax           upper pressure limit of the pressure grid [Pa]
;           prangesteps         number of pressure grid points
;           prangefile          input pressure gride file 
;
; Calls:    FUNCTION  addfields2vector
;           FUNCTION  ckeck_arts_tag
;           FUNCTION  write_arts_control_file_tag
;           FUNCTION  write_arts_control_file_contdescription
;           FUNCTION  write_arts_control_file_lineshape
;           FUNCTION  write_arts_control_file_vmr
;           PROCEDURE list_arts_tag
;           FUNCTIONS/PROCEDURES from aii_checks.pro
;
; Output:   flag                0=ok, the arts control file is correctly created
;                               1=error, the arts control file is not correctly created
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
PRO CreateArtsControlFile, flag=flag, $
                           debug=debug, $
                           artsjobpath=artsjobpath, $
                           artspath=artspath, $
                           controlfile=controlfile, $
                           tags=tags, $
                           tag_models=tag_models, $
                           tag_userparameters=tag_userparameters, $
                           catname=catname, $
                           catformat=catformat, $
                           catfmin=catfmin, $
                           catfmax=catfmax, $
                           lineshapes=lineshapes,$
                           frangemin=frangemin, $
                           frangemax=frangemax, $
                           frangesteps=frangesteps, $
                           frangefile=frangefile, $
                           ptzfile=ptzfile, $
                           vmrtagnames=vmrtagnames, $
                           vmrfilenames=vmrfilenames,$ 
                           vmrbasename=vmrbasename, $
                           vmrscenarioname=vmrscenarioname, $
                           prangemin=prangemin, $
                           prangemax=prangemax, $
                           prangesteps=prangesteps, $
                           prangefile=prangefile
;
; ---- SET DEGREE OF COMMENTS -----------------------------
IF NOT KEYWORD_SET(debug) THEN debug = 0 
IF ((debug NE 0) AND $
    (debug NE 1) AND $
    (debug NE 2) AND $
    (debug NE 3)) THEN BEGIN
    debug = 2
    print, ' CreateArtsControlFile> degree of debugging: 0=no messages'
    print, ' CreateArtsControlFile>                      1=only few important messages'
    print, ' CreateArtsControlFile>                      2=normal important messages'
    print, ' CreateArtsControlFile>                      3=detailed important messages'
    print, ' CreateArtsControlFile> Warning, you selected a debug level which is not supported!'
    print, ' CreateArtsControlFile> Set the debug level to debug=',debug
ENDIF
;
; ---- OUTPUT FLAG ----------------------------------------
flag = 1 ; 0=ok, 1=error
;
; ---- ARTS JOB SUBDIRECTORY ------------------------------
IF NOT KEYWORD_SET(artsjobpath) THEN BEGIN 
    spawn, 'pwd', artsjobpath
    artsjobpath = artsjobpath+'/'
ENDIF
;
; ---- OUTPUT FILE ----------------------------------------
IF (KEYWORD_SET(controlfile)) THEN BEGIN 
    ffile =  artsjobpath+controlfile
ENDIF ELSE BEGIN
    ffile = artsjobpath+'AutoCreateArtsControlFile.arts'
ENDELSE
IF (debug GT 0) THEN print, ' CreateArtsControlFile> arts control file with name: >>'+ffile+'<<' 
;
; ---- CHECKS ---------------------------------------------
filevec = FINDFILE(ffile)
IF (N_ELEMENTS(filevec) GT 1) THEN BEGIN
    IF (debug GT 1) THEN print, ' CreateArtsControlFile> delete arts control file(s): >>'+filevec[i]+'<<' 
    FOR i = 0, N_ELEMENTS(filevec)-1 DO spawn, 'rm -f ', filevec[i]
ENDIF
;
; ---- OPEN OUTPUT FILE -----------------------------------
openw, funit, ffile, ERROR=err, /get_lun
IF (err NE 0) THEN BEGIN
    print,' CreateArtsControlFile> ERROR: output file >>'+ffile+'<< can not be created!'
    print,' CreateArtsControlFile> error code for output file:',err
    print,' CreateArtsControlFile> error message:', !ERR_STRING
    print,' CreateArtsControlFile> STOP HERE!'
    FREE_LUN, funit & STOP
ENDIF
;
; ---- WRITE ARTS CONTROL FILE ----------------------------
; a) header information
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write header information in arts control file'
printf, funit, FORMAT='(A49)','#  ******** temporary arts control file *********'
printf, funit, FORMAT='(A49)','#  H2O continuum absorption parameter set fitting'
printf, funit, FORMAT='(A1)', '#'
;
; b) arts tags definition
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write arts tags definition in arts control file'
if (write_arts_control_file_tag(funit, tags) NE 0) THEN BEGIN
    print,' CreateArtsControlFile> ERROR: tag information can not be writen into output file!'
    FREE_LUN, funit & STOP
ENDIF
;
; c) special tag information
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write special tag information in arts control file'
IF (write_arts_control_file_contdescription(funit, tags, $
                                            tag_models, $
                                            tag_userparameters) NE 0) THEN BEGIN
    print,' CreateArtsControlFile> ERROR: tag information can not be writen into output file!'
    FREE_LUN, funit & STOP
ENDIF
;
; d) line catalog information
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write line catalog information in arts control file'
IF (KEYWORD_SET(catname)) THEN BEGIN
    IF (N_ELEMENTS(check_input_file(catname, '')) NE 1) THEN BEGIN
        print,' CreateArtsControlFile> ERROR: line cataloge file >>'+catname+'<< not found!'
        print,' CreateArtsControlFile> catname=',catname
        FREE_LUN, funit & STOP
    ENDIF ELSE BEGIN
        printf, funit, 'lines_per_tgReadFromCatalogues{'
        printf, funit, '  filenames = [ "'+catname+'" ]'
        printf, funit, '  formats   = [ "'+catformat+'" ]'
        printf, funit, '  fmin      = [ '+string(catfmin)+' ]'
        printf, funit, '  fmax      = [ '+string(catfmax)+' ]'
        printf, funit, '}'
    ENDELSE
ENDIF ELSE BEGIN
    printf, funit, 'lines_per_tgSetEmpty{}'
ENDELSE
printf, funit, FORMAT='(A1)', '#'
;
; e) line shape information for each tag
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write line shape information for each tag in arts control file'
IF (write_arts_control_file_lineshape(funit,$
                                      tags, $
                                      lineshapes) NE 0) THEN BEGIN
    print,' CreateArtsControlFile> ERROR: line shape information not correct!'
    print,lineshapes
    FREE_LUN, funit & STOP
ENDIF
printf, funit, FORMAT='(A1)', '#'
;
; f) frequency range of calculation
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write frequency range of calculation in arts control file'
IF KEYWORD_SET(frangefile) THEN BEGIN
    IF (KEYWORD_SET(frangemin) AND $
        KEYWORD_SET(frangemax) AND $
        KEYWORD_SET(frangesteps)) THEN BEGIN
        print,' CreateArtsControlFile> ERROR: ambiguity in frequency grid information!'
        print,' CreateArtsControlFile> You have specified the variables "frangefile"'
        print,' CreateArtsControlFile> and in addition one or more of the variables'
        print,' CreateArtsControlFile> "frangemin", "frangeamx", "frangesteps".'
        print,' CreateArtsControlFile> BUT you should use eider "frangefile" OR'
        print,' CreateArtsControlFile> "frangemin", "frangeamx", "frangesteps"!'
        FREE_LUN, funit & STOP
    ENDIF
    IF (N_ELEMENTS(check_input_file(frangefile, '')) EQ 1) THEN BEGIN
      printf, funit, 'VectorReadAscii(f_mono){"'+frangefile+'"}'
    ENDIF ELSE BEGIN
        print,' CreateArtsControlFile> ERROR: frequency grid information not correct!'
        print,' CreateArtsControlFile> You have not specified the variable "frangefile"'
        print,' CreateArtsControlFile> in a correct way. File >>',frangefile,'<< not found'
        print,' CreateArtsControlFile> (frangefile=path/file)'
        FREE_LUN, funit & STOP
    ENDELSE 
ENDIF ELSE BEGIN
    IF (KEYWORD_SET(frangemin) AND $
        KEYWORD_SET(frangemax) AND $
        KEYWORD_SET(frangesteps)) THEN BEGIN
        printf, funit, 'VectorNLinSpace(f_mono){'
        printf, funit, '        start =    '+string(frangemin)
        printf, funit, '        stop  =    '+string(frangemax)
        printf, funit, '        n     =    '+string(frangesteps)
        printf, funit, '}'
    ENDIF ELSE BEGIN
        print,' CreateArtsControlFile> ERROR: frequency grid information not correct!'
        print,' CreateArtsControlFile> You have not specified the variable "frangefile" or'
        print,' CreateArtsControlFile> the variables "frangemin", "frangeamx", and "frangesteps".'
        print,' CreateArtsControlFile> "frangefile":',frangefile
        print,' CreateArtsControlFile> "frangemin","frangemax","frangestep":',$
              frangemin,frangemax,frangestep
        FREE_LUN, funit & STOP
    ENDELSE
ENDELSE
printf, funit, FORMAT='(A1)', '#'
;
; g) Pressure[Pa] Temperature[K] Altitude[m] input file
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write P-T-alt file information in arts control file'
IF (KEYWORD_SET(ptzfile)) THEN BEGIN
    IF (N_ELEMENTS(check_input_file(ptzfile, '')) EQ 1) THEN BEGIN
        printf, funit, 'MatrixReadAscii (raw_ptz)'
        printf, funit, '        {"'+ptzfile+'"}'
;       /pool/lookup2/arts-data/atmosphere/fascod/midlatitude-summer.tz.aa
    ENDIF ELSE BEGIN
        print,' CreateArtsControlFile> ERROR: P-T-alt file information not correct!'
        print,' CreateArtsControlFile> There are ',check_input_file(ptzfile, ''),' files'
        print,' CreateArtsControlFile> with the specified input file name for'
        print,' CreateArtsControlFile> variable "ptzfile". STOP!'
        FREE_LUN, funit & STOP
    ENDELSE
ENDIF ELSE BEGIN
    print,' CreateArtsControlFile> ERROR: P-T-alt file information not correct!'
    print,' CreateArtsControlFile> You have not specified the variable "ptzfile". STOP!'
    FREE_LUN, funit & STOP
ENDELSE
printf, funit, FORMAT='(A1)', '#'
;
; h) VMR profiles taken from the following input files
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write VMR profile file information in arts control file'
IF (KEYWORD_SET(vmrscenarioname)) THEN BEGIN
    IF (KEYWORD_SET(vmrtagnames) OR KEYWORD_SET(vmrfilenames) OR $
        KEYWORD_SET(vmrbasename)) THEN BEGIN
        print,' CreateArtsControlFile> ERROR: ambiguity of VMR source files!'
        print,' CreateArtsControlFile> You have to provide the scenario with' 
        print,' CreateArtsControlFile> the variable "vmrscenarioname"'
        print,' CreateArtsControlFile> (for the arts method "raw_vmrsReadFromScenario")'
        print,' CreateArtsControlFile> OR'
        print,' CreateArtsControlFile> the three variables "vmrtagnames", "vmrfilenames", "vmrbasenam"'
        print,' CreateArtsControlFile> (for arts method "raw_vmrsReadFromFiles")'
        print,' CreateArtsControlFile> BUT NOT all four variables!'
        FREE_LUN, funit & STOP
    ENDIF
    printf, funit, 'raw_vmrsReadFromScenario'
    printf, funit, '         {"'+vmrscenarioname+'"}'
ENDIF ELSE BEGIN
    IF (KEYWORD_SET(vmrtagnames)) THEN BEGIN
        IF (write_arts_control_file_vmr(funit, tags, vmrtagnames, $
                                        vmrfilenames, vmrbasename) NE 0) THEN BEGIN
            print,' CreateArtsControlFile> ERROR: VMR profile file information not correct!'
            print,' CreateArtsControlFile> Check the variables:'
            print,' CreateArtsControlFile> "tags", "vmrtagnames", "vmrfilenames", "vmrbasename".'
            print,' CreateArtsControlFile> variable "ptzfile". STOP!'
            FREE_LUN, funit & STOP
        ENDIF
    ENDIF
ENDELSE
printf, funit, FORMAT='(A1)', '#'
;
; i) pressure grid of calculation
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write pressure grid information in arts control file'
IF KEYWORD_SET(prangefile) THEN BEGIN
    IF (KEYWORD_SET(prangemin) OR $
        KEYWORD_SET(prangemax) OR $
        KEYWORD_SET(prangesteps)) THEN BEGIN
        print,' CreateArtsControlFile> ERROR: ambiguity in pressure grid information!'
        print,' CreateArtsControlFile> You have specified the variables prangefile'
        print,' CreateArtsControlFile> and in addition one or more of the variables'
        print,' CreateArtsControlFile> "prangemin", "prangeamx", "prangesteps".'
        print,' CreateArtsControlFile> BUT you should use eider prangefile OR'
        print,' CreateArtsControlFile> "prangemin", "prangeamx", "prangesteps"!'
        FREE_LUN, funit & STOP
    ENDIF
    IF (N_ELEMENTS(check_input_file(prangefile, '')) EQ 1) THEN BEGIN
      printf, funit, 'VectorReadAscii(p_abs){"'+prangefile+'"}'
    ENDIF ELSE BEGIN
        print,' CreateArtsControlFile> ERROR: pressure grid information not correct!'
        print,' CreateArtsControlFile> You have not specified the variable "prangefile"'
        print,' CreateArtsControlFile> in a correct way. File >>',prangefile,'<< not found'
        print,' CreateArtsControlFile> (prangefile=path/file)'
        FREE_LUN, funit & STOP
    ENDELSE
ENDIF ELSE BEGIN
    IF (KEYWORD_SET(prangemin) AND $
        KEYWORD_SET(prangemax) AND $
        KEYWORD_SET(prangesteps)) THEN BEGIN
        printf, funit, 'VectorNLogSpace(p_abs){ '
        printf, funit, '        start =    '+string(prangemin, FORMAT='(E12.6)')
        printf, funit, '        stop  =    '+string(prangemax, FORMAT='(E12.6)')
        printf, funit, '        n     =    '+string(prangesteps)
        printf, funit, '}'
    ENDIF ELSE BEGIN
        print,' CreateArtsControlFile> ERROR: pressure grid information not correct!'
        print,' CreateArtsControlFile> You have not specified the variable prangefile or'
        print,' CreateArtsControlFile> the variables "prangemin", "prangeamx", and "prangesteps".'
        print,' CreateArtsControlFile> "prangefile":',prangefile
        print,' CreateArtsControlFile> "prangemin","prangemax","prangestep":',$
              prangemin,prangemax,prangestep
        FREE_LUN, funit & STOP
    ENDELSE
ENDELSE
printf, funit, FORMAT='(A1)', '#'
;
; j) read the atmospheric condition into arts
IF (debug GT 2) THEN print, ' CreateArtsControlFile> write arts method "AtmFromRaw{}" in arts control file'
printf, funit, 'AtmFromRaw{}' 
printf, funit, '#'
;
; k) calculate the absorption
; additionally the profiles of H2O and N2 are needed seperately
; in any case of the considered absorption calculation
printf, funit, 'h2o_absSet{}'
printf, funit, '#'
printf, funit, 'n2_absSet{}'
printf, funit, '#'
printf, funit, 'absCalc{}'
;printf, funit, 'ArrayOfMatrixWriteAscii (abs_per_tg) {""}'
printf, funit, 'MatrixWriteAscii (abs) {""}' 
;
; ---- SET FLAG TO OK -------------------------------------
flag = 0
;
; ---- END ------------------------------------------------
create_arts_control_file_ende:
FREE_LUN, funit
;
END
;
; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
