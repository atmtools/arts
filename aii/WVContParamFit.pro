; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
; ########################### INTERNAL FUNCTIONS ########################### 
;
;
;
; **************************************************************************
; Name:     WVCONTABSCALC
;
; Purpose:  calculate water vapor continuum absorption according to
;           Rosenkranz, Radio Science, 1998, parameterization.
;
; Inputs:   scalar  PWV    H2O pressure             [hPa]
;           scalar  PF     foreign pressure         [hPa]
;           scalar  F      frequency                [GHz]
;           scalar  THETA  T_ref/T                  [1]
;           scalar  CS     self cont. coeff.        [dB/km/hPa2/GHz2]
;           scalar  XS     self cont. T exponent    [1]
;           scalar  CF     foreign cont. coeff.     [dB/km/hPa2/GHz2]
;           scalar  XF     foreign cont. T exponent [1]
;
; Output:   scalar         continuum absorption [Np/m] according to
;                          Rosenkranz parameterization, Radio Science, 19998
;
; History:  2002-08-23    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION WVCONTABSCALC, PWV, PF, F, THETA, CS, XS, CF, XF 
;
; absorption in dB/km
ABS = (F * F) * (THETA * THETA * THETA) *      $
      ( (CS * THETA^XS * PWV * PWV) +          $
        (CF * THETA^XF * PWV * PF) )
;
; absorption in Np/m
ABS = 2.3025851E-04 * ABS
;
RETURN, ABS
END
;
;
; **************************************************************************
;
;
FUNCTION GetDatum
;
;; get datum to write it on the top of the plot:
spawn,'date +"%y"',year
spawn,'date +"%m"',month
spawn,'date +"%d"',day
spawn,'date +"%H"',hour
spawn,'date +"%M"',minute
spawn,'date +"%S"',second
spawn,'whoami',who
datum = who+':20'+string(year, FORMAT='(A2)')+'-'+$
string(month, FORMAT='(A2)')+'-'+$
string(day, FORMAT='(A2)')+'/'+$
string(hour, FORMAT='(A2)')+':'+$
string(minute, FORMAT='(A2)')+':'+$
string(second, FORMAT='(A2)')
RETURN, datum
END
;
;
; ********************************************************************
; Student t-distribution for double sided test of 95% 
; significance if both sided and 97.5 for single sided
; Source:
; http://www.itl.nist.gov/div898/handbook/mpc/section3/mpc3652.htm
;
FUNCTION ttest095, i

if ((i LT 1) OR (i GT 89)) then begin
    print,'ttest095> ERROR!  number of freedom is out of range!'
    print,'(1<n<=100) n =',i
    return, -999.99
endif


tdist = [ 0.000,   12.706,   4.303,   3.182,   2.776,   $
          2.571,    2.447,   2.365,   2.306,   2.262,   $
          2.228,    2.201,   2.179,   2.160,   2.145,   $
          2.131,    2.120,   2.110,   2.101,   2.093,   $
          2.086,    2.080,   2.074,   2.069,   2.064,   $
          2.060,    2.056,   2.052,   2.048,   2.045,   $
          2.042,    2.040,   2.037,   2.035,   2.032,   $
          2.030,    2.028,   2.026,   2.024,   2.023,   $
          2.021,    2.020,   2.018,   2.017,   2.015,   $
          2.014,    2.013,   2.012,   2.011,   2.010,   $
          2.009,    2.008,   2.007,   2.006,   2.005,   $
          2.004,    2.003,   2.002,   2.002,   2.001,   $
          2.000,    2.000,   1.999,   1.998,   1.998,   $
          1.997,    1.997,   1.996,   1.995,   1.995,   $
          1.994,    1.994,   1.993,   1.993,   1.993,   $
          1.992,    1.992,   1.991,   1.991,   1.990,   $
          1.990,    1.990,   1.989,   1.989,   1.989,   $
          1.988,    1.988,   1.988,   1.987,   1.987]

;,   $
;          1.987,    1.986,   1.986,   1.986,   1.986,   $
;          1.985,    1.985,   1.985,   1.984,   1.984,   $
;          1.984]

RETURN, tdist[i]
END
;
;
; **************************************************************************
; Name:     MYLINREGRESS
;
; Purpose:  performs linear regression of input vectors x andy
;
; Inputs:   vector        X        independent variable
;                         Y        dependent variable
;                         Weights  wheight of y
;
; Output:   vector        parameters a and b of y = a*x+b        
;                         additional output as keyowrds
;
; History:  2002-08-15    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION MYLINREGRESS, X, Y,              $
                       weights=weights,   $ ; weight of y
                       yfit=yfit,         $ ; yfit_i = a*x_i + b
                       mean2sep=mean2sep, $ ; mean quad. separation
                       confa=confa,       $ ; 95% confidence intervall of a
                       confy=confy,       $ 
                       hypoa=hypoa

if (N_ELEMENTS(X) NE N_ELEMENTS(Y)) then begin
    print,'MYLINREGRESS> ERROR!  wrong size of input vectors X, Y'
    print,'  N_ELEMENTS(X):',N_ELEMENTS(X),'N_ELEMENTS(Y):',N_ELEMENTS(Y)
    print,'  RETURN without calculation!'
    return, [0,0]
endif

;; number of elements in inpit vectors
N = N_ELEMENTS(X)

;; independent variable vector
xmom  = MOMENT(X, /DOUBLE)
xmean = xmom[0]
xvar  = xmom[1]

;; dependent variable vector
wn    = dblarr(N)
if keyword_set(weights) then begin
    sw = TOTAL(weights)
    sy = 0.000
    for i = 0,N-1 do wn[i] = weights[i] / sw
    for i = 0,N-1 do sy = sy + (wn[i]*Y[i])
    ymean = sy / DOUBLE(N)
    svy   = 0.000
    for i = 0,N-1 do svy = svy + ((wn[i]*Y[i])-ymean)^2
    yvar  = svy / DOUBLE(N-1)
endif else begin
    ymom  = MOMENT(Y, /DOUBLE)
    ymean = ymom[0]
    yvar  = ymom[1]
    wn[*] = 1.0000E0 / DOUBLE(N)
endelse


SXY = 0.000E0
for i = 0,N-1 do begin
    SXY = SXY + ( (X[i]-xmean) * (Y[i]-ymean) )
endfor
SXY = SXY / DOUBLE(N-1)

;; linear regression parameters for straight line y = a*x + b
a = SXY / xmom[1]
b = ymean - (a * xmean)
coeff = [a, b]


;;  --------------- additional information ---------------


if keyword_set(yfit) then begin
    yfit = dblarr(N)
    for i = 0,N-1 do begin
        yfit[i] = (a * X[i]) + b
    endfor
endif 


if keyword_set(mean2sep) then begin
;;  mean square separation between data and model
;;  the deviations in y of data amd regression line
;;  are squared and summed up. Afterwards the square 
;;  root is taken and the result devided by the number 
;;  of data points.
    dy   = dblarr(N)
    sdy  = 0.00E0
    sdy2 = 0.00E0
    for i = 0,N-1 do begin
        dy[i] = ( Y[i] - ((a * X[i]) + b) )
        sdy   = sdy + dy[i]
        sdy2  = sdy2 + (dy[i]*dy[i])
    endfor
    mean2sep = SQRT(sdy2) / DOUBLE(N)
endif 


if keyword_set(confa) then begin
;; 95% confidence intervall of the slope parameter a:
;; true a withing 95% in the interval [a-confa, a+confa]
   d = (N-1) * ( ymom[1] - (a * xmom[1]) )
   c = ttest095((N-2))
   print, 'ttest095(',N-2,')=',c 
   if (c GT 0.00) then begin
       confa = c * SQRT(d) / SQRT(xmom[1]) / SQRT(DOUBLE((N-1)*(N-2)))
   endif else begin
       print, 'WARNING!  t-distribution gave wrong number!'
       confa = -999.99
   endelse
endif 



if keyword_set(confy) then begin
;; 95% confidence intervall of the y values
;; true value withing 95% in the interval [y_i-confa_i, y_i+confa_i]
    confy = dblarr(N)
    d = DOUBLE(N-1) * ( ymom[1] - (a * xmom[1]) )
    c = ttest095((N-2))
    if (c GT 0.00) then begin
        for i = 0,N-1 do begin
            h2 = 1.0000E0 / DOUBLE(N) + $
              ( X[i]-xmean )^2 / ( (DOUBLE(N)-1) * xvar )
            confy[i] = c * SQRT(h2 * d) / SQRT(DOUBLE(N-2))
        endfor
    endif else begin
        print, 'WARNING!  t-distribution gave wrong number!'
        confy[*] = -999.99
    endelse
endif 



if keyword_set(hypoa) then begin
;; test the hypothesis that the slope parameter is 
;, significance number alpha is 2.5% 
    hypoa = 1.000
    d = DOUBLE(N-1) * ( ymom[1] - (a * xmom[1]) )
    c = ttest095((N-2))
    if (c GT 0.00) then begin
        hypo = SQRT(xvar)                        * $
          SQRT((DOUBLE(N)-1)*(DOUBLE(N)-2)) * $
          (a-hypoa) / (SQRT(d))
        if (hypo LE c) then hypoa = 0.000
    endif else begin
        print, 'WARNING!  t-distribution gave wrong number!'
    endelse
endif 



;; return with the linear regression parameters 
;; slope (a) and bias (b):  y = a*x + b
RETURN, coeff
END
;
;
; **************************************************************************
; Name:     read_H2O_data_file
;
; Purpose:  Reads data from a file in ARTS data format into a matrix.
;
; Inputs:   filename      full file name
;
; Output:   matrix        the data matrix 
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION read_H2O_data_file, filename
;
; ---- SET OUTPUT FLAG ------------------------------------
flag = 1
;
; ---- CHECKS ---------------------------------------------
filevec = FINDFILE(filename)
IF (N_ELEMENTS(filevec) GT 1) THEN BEGIN
    print,' read_H2O_data_file> !!! ERROR: can not find data file'
    print,' read_H2O_data_file> !!! >>> multiple files found with name: ',filename,':',filevec,'! STOP!'
    stop
ENDIF
IF (N_ELEMENTS(filevec) LT 1) THEN BEGIN
    print,' read_H2O_data_file> !!! ERROR: can not find data file'
    print,' read_H2O_data_file> !!! >>> no files found with name: ',filename,'! STOP!'
    stop
ENDIF
;
; ---- DEFINITIONS ----------------------------------------
;
commentsym = '#' ; symbol for comment lines
datarow   = commentsym ; string which contains one line of input data
; columns entry information:
;# 0: buffer gas species                  string
dummy0 = ' '
;# 1: buffer gas pressure         [Pa]    double
dummy1 = 0.0e0
;# 2: H2O species                         string
dummy2 = ' '
;# 3: H2O pressure                [Pa]    double
dummy3 = 0.0e0
;# 4: frequency of measurement    [Hz]    double
dummy4 = 0.0e0
;# 5: termperature of measurement [K]     double
dummy5 = 0.0e0
;# 6: total measured absorption   [1/m]   double
dummy6 = 0.0e0
;# 7: data source                         integer
dummy7 = ' '
data_file_def_colums = 8 ; (0-7)
;
; ---- OPEN FILE FOR READING ------------------------------
openr, unit, filename, error=err, /GET_LUN
IF err NE 0 THEN BEGIN
    print,' read_H2O_data_file> !!! ERROR: can not open H2O absorption data file '+filename+'.'
    print,' read_H2O_data_file> !!! Please check the directory and file names! STOP!'
    print,' read_H2O_data_file> !!! Error message:', !ERR_STRING
    FREE_LUN, unit
    stop
ENDIF
;
; ---- READ HEADER INFORMATION ----------------------------
; Read until line does not begin with #
WHILE STRMID(STRTRIM(datarow, 1), 0, 1) EQ commentsym DO BEGIN 
    readf, unit, datarow
ENDWHILE
IF datarow EQ '' THEN BEGIN
    print, '!!! ERROR: on input line:',datarow
    print, '!!! >>> blank lines are not allowed! STOP!'
    FREE_LUN, unit
    stop
ENDIF
;
; ---- READ ROW/COLUMN numbers ----------------------------
nrows = 0
ncols = 0
;readf, unit, datarow
p = str_sep(datarow, ' ')
IF (N_ELEMENTS(p) NE 2) THEN BEGIN
    print,'!!! ERROR: array size: # rows=',p[0],', # colums=',p[1]
    print,'!!! >>> wrong array size? STOP!'
    FREE_LUN, unit & stop
ENDIF
reads, p[0], nrows
reads, p[1], ncols
IF ncols NE data_file_def_colums THEN BEGIN
    print,'!!! ERROR: no of columns: # colums=',p[1]
    print,'!!! >>> wrong size? STOP!'
    FREE_LUN, unit & stop
ENDIF
;
print, 'input data array: # rows=',nrows,', # colums=',ncols
;
; ---- DEFINE OUTPUT STRUCTURE ARRAY ----------------------
; BuffName:  buffer gas species
; pd:        buffer gas pressure [Pa]
; WVName:    water vapor species
; pwv:       water vapor pressure [Pa]
; f:         frequnce of measurement [Hz]
; T:         temperature [K]
; abs:       measured total absorption [1/m] 
; LineMod:   line absorption model name
; labs:      calculated line absorption [1/m]
; cabs:      calculated line absorption [1/m]
rec = {MWVDATA, BuffName:'', pd:0.0e0, WVName:'', pwv:0.0e0, f:0.0e0, T:0.0e0, abs:0.0e0, source:'', $
                LineMod:'', labs:0.0e0, cabs:0.0e0}
measurements = REPLICATE(rec, nrows)
;
; ---- READ INPUT DATA ------------------------------------
; a) Read the matrix row by row
FOR r = 0, nrows - 1 DO BEGIN
    readf, unit, datarow
    p = STR_SEP(STRCOMPRESS(datarow), ',')
    IF N_ELEMENTS(p) NE ncols THEN BEGIN
        print, '!!! ERROR: Wrong number of columns detected in data file! STOP!'
        FOR i=0,N_ELEMENTS(p)-1 DO print, 'input line vecor p[',i,']=>>'+p[i]+'<<'
        stop
    ENDIF
    reads, p[0], dummy0
    reads, p[1], dummy1
    reads, p[2], dummy2
    reads, p[3], dummy3
    reads, p[4], dummy4
    reads, p[5], dummy5
    reads, p[6], dummy6
    reads, p[7], dummy7
    q = STR_SEP(dummy0, '"')
    measurements[r].BuffName = q[1]
    measurements[r].pd       = dummy1
    q = STR_SEP(dummy2, '"')
    measurements[r].WVName   = q[1]
    measurements[r].pwv      = dummy3
    measurements[r].f        = dummy4
    measurements[r].T        = dummy5
    measurements[r].abs      = dummy6
    measurements[r].source   = dummy7
    measurements[r].labs     = 0.0e0
    measurements[r].cabs     = 0.0e0
ENDFOR
;
;
; ---- CHECK IF THE FILE IS NOW FINISHED -----------
datarow = ''
WHILE NOT eof(unit) DO BEGIN
  readf, unit, datarow
  s = strcompress(datarow, /remove_all)
  IF s NE '' THEN BEGIN
    print, '!!! ERROR: there is some additional data at the end of the file:'
    print, '"',s,'"'
    print,'!!! >>> wrong array size? STOP!'
    FREE_LUN, unit & stop
  ENDIF
ENDWHILE
;
; ---- CLOSE THE FILE ------------------------------
flag = 0
free_lun, unit  
;
; ---- RETURN MATRIX -------------------------------
read_H2O_data_file_ende:
RETURN, measurements
;
END
;
; **************************************************************************
; Name:     write_ptz_file
;
; Purpose:  Reads data from a file in ARTS data format into a matrix.
;
; Inputs:   filename      full file name
;
; Output:   matrix        the data matrix 
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION write_ptz_file, artsjobpath, controlfile, ptot, T
;
;
; ---- P-T-z FILE NAME ------------------------------------
ptzfile = artsjobpath+controlfile+'_ptz.aa'
;
; ---- CHECKS ---------------------------------------------
filevec = FINDFILE(ptzfile)
IF (N_ELEMENTS(filevec) GT 1) THEN BEGIN
    IF (debug GT 1) THEN print, ' write_ptz_file> delete arts P-T-z file: >>'+filevec[i]+'<<' 
    FOR i = 0, N_ELEMENTS(filevec)-1 DO spawn, 'rm -f ', filevec[i]
ENDIF
;
; ---- OPEN OUTPUT FILE -----------------------------------
openw, unit, ptzfile, ERROR=err, /get_lun
IF (err NE 0) THEN BEGIN
    print,' write_ptz_file> ERROR: output file >>'+ptzfile+'<< can not be created!'
    print,' write_ptz_file> error code for output file:',err
    print,' write_ptz_file> error message:', !ERR_STRING
    FREE_LUN, unit
    goto, write_ptz_file_ende
ENDIF
;
; ---- WRITE OUTPUT FILE ----------------------------------
printf, unit, '# ARRAY dimension (=1)'
printf, unit, '# MATRIX dimensions (always 3 columns)'
printf, unit, '# Pressure[Pa] Temperature[K] Altitude[m]'
printf, unit, '# This file is created by the IDL function write_ptz_file'
printf, unit, '1'
printf, unit, '2  3'
printf, unit, string(2.000*(ptot), FORMAT='(E12.6)')+' '+string(T, FORMAT='(E12.6)')+' 0.000000e+00'
printf, unit, string(0.500*(ptot), FORMAT='(E12.6)')+' '+string(T, FORMAT='(E12.6)')+' 1.000000e+00'
;
; ---- END OF FUNCTION ------------------------------------
write_ptz_file_ende:
FREE_LUN, unit
RETURN, ptzfile
END
;
; **************************************************************************
; Name:     write_vmr_file
;
; Purpose:  Reads data from a file in ARTS data format into a matrix.
;
; Inputs:   filename      full file name
;
; Output:   matrix        the data matrix 
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION write_vmr_file, artsjobpath, controlfile, tags, ptot, p
;
; ---- CHECKS ---------------------------------------------
;IF (ABS(N_ELEMENTS(tags)-N_ELEMENTS(p)) GT 0) THEN BEGIN
;    print, ' write_vmr_file> ERROR in size of tags and pressure detected!'
;    print, ' write_vmr_file> pressure and tag vectors must have the same size!'
;    print, ' write_vmr_file> tag vector size     :',N_ELEMENTS(tags)
;    print, ' write_vmr_file> pressure vector size:',N_ELEMENTS(p)
;    FREE_LUN, unit & STOP
;ENDIF
;
; ---- DELETE EXISTING FILES ------------------------------
vmrfile = strarr(N_ELEMENTS(tags))
FOR k = 0, N_ELEMENTS(tags)-1 DO BEGIN
    vmrfile[k] = artsjobpath+controlfile+'_vmr_'+tags[k]+'.aa'
    filevec = FINDFILE(vmrfile[k])
    IF (N_ELEMENTS(filevec) GT 1) THEN BEGIN
        IF (debug GT 1) THEN print, ' write_ptz_file> delete arts P-T-z file: >>'+filevec[i]+'<<' 
        FOR i = 0, N_ELEMENTS(filevec)-1 DO spawn, 'rm -f ', filevec[i]
    ENDIF
ENDFOR
;
; ---- LOOP OVER ALL TAGS -----------------------------------
FOR i = 0,  N_ELEMENTS(tags)-1 DO BEGIN
;   ---- OPEN OUTPUT FILE -----------------------------------
    openw, unit, vmrfile[i], ERROR=err, /get_lun
    IF (err NE 0) THEN BEGIN
        print,' write_ptz_file> ERROR: output file >>'+vmrfile[k]+'<< can not be created!'
        print,' write_ptz_file> error code for output file:',err
        print,' write_ptz_file> error message:', !ERR_STRING
        FREE_LUN, unit & goto, write_ptz_file_jump
    ENDIF
;
;   ---- CALCULATIONS ---------------------------------------
    vmr = p[i] / ptot
;
;   ---- WRITE OUTPUT FILE ----------------------------------
    printf, unit, '# ARRAY dimension (=1)'
    printf, unit, '# MATRIX dimensions (always 2 columns)'
    printf, unit, '# Pressure[Pa] VMR[absolute number]'
    printf, unit, '# This file is created by the IDL function write_vmr_file'
    printf, unit, '1'
    printf, unit, '2  2'
    printf, unit, string(2.000*(ptot), FORMAT='(E12.6)')+' '+string(vmr, FORMAT='(E12.6)')
    printf, unit, string(0.500*(ptot), FORMAT='(E12.6)')+' '+string(vmr, FORMAT='(E12.6)')
;
    write_ptz_file_jump:
    FREE_LUN, unit 
ENDFOR
;
; ---- END OF FUNCTION ------------------------------------
write_vmr_file_ende:
RETURN, vmrfile
END
;
; **************************************************************************
; Name:     read_lineabs_data
;
; Purpose:  Reads line absorption data from arts output file
;
; Inputs:   filename      full file name with extension
;
; Output:   labs          line absorption
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION read_lineabs_data, artsjobpath, filename
;
;
; ---- OUTPUT VARIABLE ------------------------------------
labs = 0.0e0
;
; ---- CHECKS ---------------------------------------------
file = artsjobpath+filename
filevec = FINDFILE(file)
IF (N_ELEMENTS(filevec) GT 1) THEN BEGIN
    print,' read_lineabs_data> !!! ERROR: can not find data file'
    print,' read_lineabs_data> !!! >>> multiple files found with name: ',file,':',filevec,'! STOP!'
    STOP
ENDIF
IF (N_ELEMENTS(filevec) LT 1) THEN BEGIN
    print,' read_lineabs_data> !!! ERROR: can not find data file'
    print,' read_lineabs_data> !!! >>> no files found with name: ',file,'! STOP!'
    STOP
ENDIF
;
; ---- READ FILE -------------------------------------------
print, ' read_lineabs_data> filename:',file
absmat = aa_read_general(file)
;
; ---- EXTRACT LINE ABSORPTION FORM MATRIX -----------------
labs = absmat[0,0]
;
; ---- END OF FUNCTION ------------------------------------
read_lineabs_data_ende:
RETURN, labs
END
;
; **************************************************************************
; Name:     pwrcontiform1
;
; Purpose:  needed for the fit of data with fit procedure svdfit
;
; Inputs:   X
;           M
;
; Output:   derivatives
;
; History:  2002-01-16    Thomas Kuhn, iup Bremen
;
; **************************************************************************
;
FUNCTION pwrcontiform1, X ,M
;
; y = a0 + a1*x
; 
;      [ df/da0,   df/da1]
return,[ [1.0e0],    [X] ]
;
END
;
; **************************************************************************
; Name:     pwrcontiform2
;
; Purpose:  needed for the fit of data with fit procedure lmfit
;
; Inputs:   X
;           M
;
; Output:   derivatives
;
; History:  2002-01-16    Thomas Kuhn, iup Bremen
;
;
; **************************************************************************
;
FUNCTION pwrcontiform2, X ,A
;
; y = a0 + a1*x
;
;      [ f(x),          df/da0,   df/da1]
RETURN,[ [A[0]+A[1]*X], [1.0E0] , [X] ]
;
END
;
;--------------------------------------------------------------------------
;
; **************************************************************************
; Name:     pwrcontiform3
;
; Purpose:  needed for the fit of data with fit procedure lmfit
;
; Inputs:   X
;           M
;
; Output:   derivatives
;
; History:  2002-01-16    Thomas Kuhn, iup Bremen
;
;
; **************************************************************************
;
FUNCTION pwrcontiform3, X ,A
;---------------------------
;
; y = a0 * x^a1
;
;      [ f(x),          df/da0,    df/da1          ]
RETURN,[ [A[0]*X^A[1]], [X^A[1]] , [A[0]*X^A[1]*X] ]
;
END
;
; **************************************************************************
; Name:     find_entry_in_vec
;
; Purpose:  searches in a vector a specific entry
;
; Inputs:   value         Reference entry value
;           vector        Vector in which the reference value is searched
;
; Output:   index         If <value> is found in <vector>, index
;                         gives the first vector element index in
;                         which <value> was found. 
;                         The maximum difference between <value> and
;                         the vector entry is set to 0.1 per cent (=0.001)
;                         If <value> is not found in <vector>, the return
;                         value is -1
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION find_entry_in_vec, value, vector
;
; ---- SET OUTPUT FLAG ------------------------------------
i = -1
;
; ---- CHECKS ---------------------------------------------
; if the input vector is not correctly in size, give -1
IF (N_ELEMENTS(vector) LT 1) THEN BEGIN
    RETURN, i
ENDIF
;
; ---- CHECKS ---------------------------------------------
;
FOR k = 0,N_ELEMENTS(vector)-1 DO BEGIN
    IF ( ABS((value-vector[k]) / value) LT 0.001) THEN RETURN, k
ENDFOR
    
RETURN, i
END
;
;
; **************************************************************************
; Name:     get_userdate_info
;
; Purpose:  searches in a vector a specific entry
;
; Inputs:   
;
; Output:   datum         string of user name and actual date.
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION get_userdate_info, flag
;
IF (flag EQ 0) THEN BEGIN
; get datum to write it on the top of the plot:
spawn,'date +"%y"',year
spawn,'date +"%m"',month
spawn,'date +"%d"',day
spawn,'date +"%H"',hour
spawn,'date +"%M"',minute
spawn,'whoami',who
datum = who+':'+string(year, FORMAT='(A2)')+'-'+$
  string(month, FORMAT='(A2)')+'-'+string(day, FORMAT='(A2)')+'/'+$
  string(hour, FORMAT='(A2)')+':'+string(minute, FORMAT='(A2)')
ENDIF
;
spawn,'date +"%y"',year
spawn,'date +"%m"',month
spawn,'date +"%d"',day
spawn,'whoami',who
datum = who+'/'+string(year, FORMAT='(A2)')+'-'+$
  string(month, FORMAT='(A2)')+'-'+string(day, FORMAT='(A2)')
;
RETURN, datum
END
;
; **************************************************************************
; NAME:     aii_csf_plot_ps_definitions
;
; PURPOSE:  set all the relevant grahical device definitions for
;           Postscript files
;
; INPUTS:   -
;
; OUTPUT:   colors        plot color array (30 different color)
;           lstyle        line style in plots (5 different ones)
;           pstyle        symbol style in plots (13 different ones)
;           plotpos       position of the plot (normalized coordinates)
;           thick         line thickness
;
; HISTORY:  2002-01-16    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
PRO aii_csf_plot_ps_definitions, colors=colors, $
                                 lstyle=lstyle, $
                                 pstyle=pstyle, $
                                 plotpos=plotpos, $
                                 thick=thick

; --- set frame of the plot in lin-lin scalo of log-log variables ----------
; make 1 plot per page
!P.multi    = [0,1,1]
; settings for the plot
!X.MARGIN   = [0,0]
!Y.MARGIN   = [0,0]
!X.OMARGIN  = [0,0]
!Y.OMARGIN  = [0,0]
!P.REGION   = [0.0, 0.0, 0.0, 0.0]
!P.POSITION = [0.25, 0.20, 0.80, 0.9]
!P.CHARTHICK = 5.0
!P.FONT      = 1
!P.CHARSIZE  = 1.5
!X.CHARSIZE  = 1
!Y.CHARSIZE  = 1
!P.THICK = 5
!X.THICK = 5
!Y.THICK = 5
;
; --- OUTPUT VARIABLES ----------------------------------------------------
plotpos = !P.POSITION
thick   = !P.THICK
for i=0,N_ELEMENTS(colors)-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in aii_color_table
    lstyle[i] = i mod  5    ;; we have  5 different line styles
    pstyle[i] = i mod 13    ;; we have 13 different plot symbols
endfor

;
END
;
; **************************************************************************
;$Id: WVContParamFit.pro,v 1.1.2.3 2003/09/08 13:11:36 cmels Exp $
;
; Copyright (c) 1994-1998, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
;+
;NAME:
;       LINFIT
;
; PURPOSE:
;       This function fits the paired data {X(i), Y(i)} to the linear model,
;       y = A + Bx, by minimizing the chi-square error statistic. The result
;       is a two-element vector containing the model parameters,[A,B].
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = LINFIT(X, Y)
;
; INPUTS:
;       X:    An n-element vector of type integer, float or double.
;
;       Y:    An n-element vector of type integer, float or double.
;
; KEYWORD PARAMETERS:
;   CHISQ:    Use this keyword to specify a named variable which returns the
;             chi-square error statistic as the sum of squared errors between
;             Y(i) and A + BX(i). If individual standard deviations are
;             supplied, then the chi-square error statistic is computed as
;             the sum of squared errors divided by the standard deviations.
;
;  DOUBLE:    If set to a non-zero value, computations are done in double
;             precision arithmetic.
;
;    PROB:    Use this keyword to specify a named variable which returns the
;             probability that the computed fit would have a value of CHISQR
;             or greater. If PROB is greater than 0.1, the model parameters
;             are "believable". If PROB is less than 0.1, the accuracy of the
;             model parameters is questionable.
;
;    SDEV:    An n-element vector of type integer, float or double that
;             specifies the individual standard deviations for {X(i), Y(i)}.
;
;   SIGMA:    Use this keyword to specify a named variable which returns a
;             two-element vector of probable uncertainties for the model par-
;             ameters, [SIG_A,SIG_B].
;
;
; EXAMPLE:
;       Define two n-element vectors of paired data.
;         x = [-3.20, 4.49, -1.66, 0.64, -2.43, -0.89, -0.12, 1.41, $
;               2.95, 2.18,  3.72, 5.26]
;         y = [-7.14, -1.30, -4.26, -1.90, -6.19, -3.98, -2.87, -1.66, $
;              -0.78, -2.61,  0.31,  1.74]
;       Define a vector of standard deviations with a constant value of 0.85
;         sdev = replicate(0.85, n_elements(x))
;       Compute the model parameters, A and B.
;         result = linfit(x, y, chisq = chisq, prob = prob, sdev = sdev)
;       The result should be the two-element vector:
;         [-3.44596, 0.867329]
;       The keyword parameters should be returned as:
;         chisq = 11.4998, prob = 0.319925
;
; REFERENCE:
;       Numerical Recipes, The Art of Scientific Computing (Second Edition)
;       Cambridge University Press
;       ISBN 0-521-43108-5
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, September 1994
;                    LINFIT is based on the routines: fit.c, gammq.c, gser.c,
;                    and gcf.c described in section 15.2 of Numerical Recipes,
;                    The Art of Scientific Computing (Second Edition), and is
;                    used by permission.
;         Modified:  SVP, RSI, June 1996
;                    Changed SIG_AB to SIGMA to be consistant with the other
;                    fitting functions. Changed CHISQR to CHISQ in the docs
;                    for the same reason. Note that the chisqr and the SIG_AB
;                    keywords are left for backwards compatibility.
;         Modified:  GGS, RSI, October 1996
;                    Modified keyword checking and use of double precision.
;                    Added DOUBLE keyword.
;
;         2002-06-19
;         TKS copy of LINFIT of IDL resource to modidy it for fit purposes
;-
; --------------------------------------------------------------------------
; 
FUNCTION MyLinFit, x, y, chisqr = chisqr, Double = Double, prob = prob, $
                   sdev = sdev, sig_ab = sig_ab, sigma = sigma,         $
                   covab = covab, rab = rab


  ON_ERROR, 2
 
  TypeX = SIZE(X)
  TypeY = SIZE(Y)
  nX = TypeX[TypeX[0]+2]
  nY = TypeY[TypeY[0]+2]
 
  if nX ne nY then $
    MESSAGE, "X and Y must be vectors of equal length."
 
  ;If the DOUBLE keyword is not set then the internal precision and
  ;result are identical to the type of input.
  if N_ELEMENTS(Double) eq 0 then $
    Double = (TypeX[TypeX[0]+1] eq 5 or TypeY[TypeY[0]+1] eq 5)
 
  nsdev = n_elements(sdev)

 
  if nsdev eq nX then begin ;Standard deviations are supplied.
    wt = 1.0 / sdev^2
    ss = TOTAL(wt, Double = Double)                ; (15.2.4)
    sx = TOTAL(wt * x, Double = Double)            ; (15.2.4)
    sy = TOTAL(wt * y, Double = Double)            ; (15.2.4)
    t =  (x - sx/ss) / sdev                        ; (15.2.15)
    st2 = TOTAL(t^2, Double = Double)              ; (15.2.16)
    b = TOTAL(t * y / sdev, Double = Double)
  endif else if nsdev eq 0 then begin
    ss = nX + 0.0                                  ; (15.2.4)
    sx = TOTAL(x, Double = Double)                 ; (15.2.4)
    sy = TOTAL(y, Double = Double)                 ; (15.2.4)
    t = x - sx/ss                                  ; (15.2.15)
    st2 = TOTAL(t^2, Double = Double)              ; (15.2.16)
    b = TOTAL(t * y, Double = Double)
  endif else $
    MESSAGE, "sdev and x must be vectors of equal length."
 
  if Double eq 0 then begin
    st2 = FLOAT(st2) & b = FLOAT(b)                ; (15.2.16)
  endif

 
  b = b / st2                                      ; (15.2.17)
  a = (sy - sx * b) / ss                           ; (15.2.18)
  sdeva = SQRT((1.0 + sx * sx / (ss * st2)) / ss)  ; (15.2.19)
  sdevb = SQRT(1.0 / st2)                          ; (15.2.20)
 
; introduced by TKS
  covab = -sx / (ss * st2)                         ; (15.2.21)
  rab   = covab / (sdeva * sdevb)                  ; (15.2.22)

  if nsdev ne 0 then begin; if user has specified the error on y
    chisqr = TOTAL( ((y - a - b * x) / sdev)^2, Double = Double )
    if Double eq 0 then chisqr = FLOAT(chisqr)
    prob = 1 - IGAMMA(0.5*(nX-2), 0.5*chisqr)
;
  endif else begin; if user has not specified the error on y
    chisqr = TOTAL( (y - a - b * x)^2, Double = Double )
    if Double eq 0 then chisqr = FLOAT(chisqr)
    prob = chisqr * 0 + 1 ;Make prob same type as chisqr.
    sdevdat = SQRT(chisqr / (nX-2))                ; (15.1.6)
    sdeva = sdeva * sdevdat                        ; (15.1.6)
    sdevb = sdevb * sdevdat                        ; (15.1.6)
  endelse
 
  sig_ab = [sdeva, sdevb]
  sigma  = sig_ab
 
  if Double eq 0 then RETURN, FLOAT([a, b]) else RETURN, [a, b]
 
END
;
;
; **************************************************************************
; NAME:     do_PWR98_h2oh2o_cont_fit
;
; PURPOSE:  Perform a H2O continuum parameter fit according to the water
;           vapor continuum parameterization of P. W. Rosenkranz,
;           Radio Science, 33(4),  919, 1998, Radio Science, 34(4), 1025, 1999
;           The only data considered is pure H2O.
;
; INPUTS:   filename      idl structure with the all data
;
; OUTPUT:   flag          0: ok  / else: false
;           file          the file is named after the core file name
;                         of this job with the extension '.param'
;
; EXTERNAL: TeXtoIDL package (http://physweb.mnstate.edu/mcraig/TextoIDL/v2-0/default.htm)
;           CreateArtsControlFile.pro (arts/aii)
;           aii_checks.pro            (arts/aii)
;           aii_plot_file.pro         (arts/aii)
;           aii_plotsymbols.pro       (arts/aii)
;
; HISTORY:  2002-01-16    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION do_PWR98_h2oh2o_cont_fit, data, artsjobpath, controlfilenamecore, $
                                   linetotmaxratio,                        $
                                   FFITMIN, FFITMAX,                       $
                                   plotfilename, plotfileformat,           $
                                   titletext

;
;
; ---- SAVE SETTINGS -------------------------------------------------------
P_ini = !P
;
; ---- GENERAL CONSTANTS ---------------------------------------------------
COMMON UNITCONVERSION, TENLOG10_EULER, dBkm2Npm, Npm2dBkm, Hz2GHz, Pa2hPa
;
Tref = 300.0E0                    ; [K] for Theta calculation
Rlinetotmax          = 1.000      ; ratio of line to total absorption
Rlinetotmin          = 0.001      ; ratio of line to total absorption
datanlim             = 3          ; minimum of data points for fit
sigma_ratio_abs_tot  = 1.000e-1   ; Gaussian error propagation: sigma_abstot = 10.0% abs_tot 
sigma_ratio_abs_line = 5.000e-2   ; Gaussian error propagation: sigma_absl   =  5.0% abs_l
sigma_ratio_f        = 1.000e-5   ; Gaussian error propagation: sigma_f      =  0.001% f
sigma_ratio_pwv      = 1.000e-3   ; Gaussian error propagation: sigma_PH2O   =  0.1% P_H2O
sigma_ratio_T        = 5.000e-3   ; Gaussian error propagation. sigma_T      =  0.5% T
;sigma_ratio_abs_tot  = 2.000e-1  ; Gaussian error propagation: sigma_abstot = 20.0% abs_tot 
;sigma_ratio_abs_line = 1.000e-1  ; Gaussian error propagation: sigma_absl   = 10.0% abs_l
;sigma_ratio_f        = 1.000e-3  ; Gaussian error propagation: sigma_f      =  0.1% f
;sigma_ratio_pwv      = 1.000e-3  ; Gaussian error propagation: sigma_PH2O   =  0.1% P_H2O
;sigma_ratio_T        = 5.000e-3  ; Gaussian error propagation. sigma_T      =  0.5% T
IF (linetotmaxratio GE 1.00) THEN linetotmaxratio = Rlinetotmax 
;
; --- FIT VARIABLES H2O-H2O ------------------------------------------------
; these are the variables with the fit results for max. n_fits different
; fit procedures
n_fits       = 3
CS           = dblarr(3,n_fits) ; C_s for different fit routines
                                ; 0 : C_s 
                                ; 1 : x_s 
                                ; 2 : x_f (if performed) 
CS[*,*]      = 0.0e0
SIGMA_CS     = dblarr(4,n_fits) ; std. deviations for different fit routines
                                ; 0 : +sigma of C_s
                                ; 1 : -sigma of C_s
                                ; 2 : sigma of x_s
                                ; 3 : sigma of x_f (if freq. fit is performed)
SIGMA_CS[*,*]= 0.0e0
CHI2_CS      = dblarr(n_fits)   ; chi² of the different fit routines
CHI2_CS[*]   = 0.0e0
COVAB=0.000
RAB=0.000
char_fit_fun = strarr(n_fits)   ; strings, names the for different fit routines
;
; --- SORT H2O-H2O DATA ----------------------------------------------------
h2oh2o = dblarr(N_ELEMENTS(data[*].abs), 4) ; 0 : y = ln(a)
                                            ; 1 : x = ln(Theta)
                                            ; 2 : frequency
                                            ; 3 : sigma_y
j = 0
FOR i = 0, N_ELEMENTS(data[*].abs)-1 DO BEGIN
    IF ((data[i].BuffName EQ 'XX') AND (data[i].WVName EQ 'H2O'))  THEN BEGIN
;   --- check frequency range ------------------------------------------------
        IF ( ( data[i].f GE FFITMIN ) AND $
             ( data[i].f LE FFITMAX )   ) THEN BEGIN
;   --- check line to total absorption ratio ---------------------------------
        IF ( ((data[i].labs / data[i].abs) LT linetotmaxratio) AND $
             ((data[i].labs / data[i].abs) GT Rlinetotmin) ) THEN BEGIN
;         o Theta [1]
            h2oh2o[j,1] = (Tref / data[i].T)
;         o "absorption" term [dB/km]
            z = Npm2dBkm * ( data[i].abs - data[i].labs ) ; [dB/km]
;         o nominator [GHz²*hPa²]
            n = (Hz2GHz*data[i].f)^2    *  $              ; f²     [GHz²]
                (Pa2hPa*data[i].pwv)^2  *  $              ; P_H2O² [hPa²]
                (h2oh2o[j,1])^3                           ; Theta  [1]
            a = (z/n)                                     ; => a_i = (z/n)
            h2oh2o[j,0] = ALOG( a )                       ; => y_i = ln(a_i) with a_i = (z/n)
;         o calculate an error estimate for  absorption [1/m/GHz²/hPa²]
;           with Gaussian error propagation
            diffa = dblarr(5) ; differentials of a with respect to all variables
            siga  = dblarr(5) ; errors of all variables (estimated)
            diffa[0] =  1.000e0/n                              ; d a / d abs_tot
            diffa[1] = -1.000e0/n                              ; d a / d abs_l
            diffa[2] =  a * (-2.00e0/(Hz2GHz*data[i].f))       ; d a / d f 
            diffa[3] =  a * (-2.00e0/(Pa2hPa*data[i].pwv))     ; d a / d P_H2O
            diffa[4] =  a * (-3.00e0/h2oh2o[j,1])              ; d a / d Theta
            siga[0]  = sigma_ratio_abs_tot  * (Npm2dBkm*data[i].abs)
            siga[1]  = sigma_ratio_abs_line * (Npm2dBkm*data[i].labs)
            siga[2]  = sigma_ratio_f        * (Hz2GHz*data[i].f)
            siga[3]  = sigma_ratio_pwv      * (Pa2hPa*data[i].pwv)
            siga[4]  = sigma_ratio_T        * h2oh2o[j,1]
            sig_a    = 0.000e0
            FOR k = 0, N_ELEMENTS(diffa)-1 DO BEGIN
                sig_a = sig_a + ( (diffa[k])^2 * (siga[k])^2 )
;                print, k,': diffa=',diffa[k],', siga=',siga[k],', sig_a=',sig_a
            ENDFOR
            sig_a   = SQRT( sig_a )
            sig_y_p = ABS( ALOG((z/n)+sig_a) - h2oh2o[j,0] ) ; upper side
            sig_y_m = ABS( ALOG((z/n)-sig_a) - h2oh2o[j,0] ) ; lower side
            h2oh2o[j,3] = 0.500e0*(sig_y_p+sig_y_m) ; take the mean of both
;         o Ln(Theta) [1]
            h2oh2o[j,1] = -ALOG( h2oh2o[j,1] ) ; => x_i = -ln(Theta)
;         o frequency [Hz]
            h2oh2o[j,2] = Hz2GHz*data[i].f
;         o print data for control
            print,j,': x=',h2oh2o[j,1],', y=',h2oh2o[j,0],', s_y=',h2oh2o[j,3]
;         o set data counter one higher
            j = j + 1
        ENDIF
        ENDIF
    ENDIF
ENDFOR
IF (j LT datanlim) THEN BEGIN
    print,' do_PWR98_cont_fit> too few data points (#=',j,') for H2O-H2O fit.'
    print,' do_PWR98_cont_fit> jump to the end of this function!'
    goto, read_lineabs_data_ende
ENDIF
print,' do_PWR98_cont_fit> no. of H2O-H2O data points for fit=',j
h2oh2o2 = h2oh2o[0:j-1, 0:N_ELEMENTS(h2oh2o[0,*])-1] ; resize array
index   = SORT(h2oh2o2[*,1])                    ; sort in increasing ln(Theta)
vecsize = N_ELEMENTS(h2oh2o2[*,0])              ; just for simplicity
IF (N_ELEMENTS(index) NE vecsize) THEN BEGIN
    print,' do_PWR98_cont_fit> !!! ERROR: inconsistency in vector size!'
    print,' do_PWR98_cont_fit> !!! N_ELEMENTS(index):',N_ELEMENTS(index)
    print,' do_PWR98_cont_fit> !!! vecsize          :',vecsize
    goto, read_lineabs_data_ende
ENDIF
FOR i = 0, vecsize-1 DO BEGIN
    j = index[i]
    h2oh2o[i,0] = h2oh2o2[j,0] ; order the array in increasing ln(Theta)
    h2oh2o[i,1] = h2oh2o2[j,1]
    h2oh2o[i,2] = h2oh2o2[j,2]
    h2oh2o[i,3] = h2oh2o2[j,3]
    print, FORMAT='(I2,A4,F8.4,A6,F8.4,A4,F6.3,A4,F7.3)', $
           i,': y=',h2oh2o[i,0],': s_y=',h2oh2o[i,3],', x=',h2oh2o[i,1],', f=',h2oh2o[i,2]
ENDFOR
h2oh2o = h2oh2o2[0:vecsize-1, 0:N_ELEMENTS(h2oh2o[0,*])-1]      ; resize array
;
;
; ==========================={FIT H2O-H2O DATA}=============================
;
;
fitway = 1
IF (fitway EQ 1) THEN BEGIN ;  ----------- LINFIT --------------------------
;   The LINFIT function fits the paired data {xi, yi} to the 
;   linear model, y = A + Bx, by minimizing the Chi-square error statistic. 
;   The result is a two-element vector containing the model parameters [A, B]. 
; o perform fit
    sdevvec = dblarr( N_ELEMENTS(h2oh2o[0:vecsize-1,0]) )
    for isdev = 0, N_ELEMENTS(sdevvec)-1 do  $
      sdevvec[isdev] = 1.000 / ABS( h2oh2o[isdev,0] )
    COVAB = 0.000
    RAB   = 0.000
    coeff = MYLINFIT(h2oh2o[0:vecsize-1,1],    $
                   h2oh2o[0:vecsize-1,0],      $
;                   SDEV=sdevvec,               $
                   PROB=probability , $
                   SIGMA=error,       $
                   CHISQ=CHISQ2,      $
                   COVAB=COVAB,       $
                   RAB=RAB,           $
                   /DOUBLE)
    CS[0,0]       = EXP(coeff[0])                   ; fit values of C_s
    CS[1,0]       = coeff[1]                        ; fit values of x_s
    CHI2_CS[0]    = CHISQ2                          ; chi^2 of the fit
    SIGMA_CS[0,0] = CS[0,0]*(EXP( error[0])-1.0e0)  ; sigma_plus of C_s
    SIGMA_CS[1,0] = CS[0,0]*(EXP(-error[0])-1.0e0)  ; sigma_minus of C_s
    SIGMA_CS[2,0] = error[1]                        ; error in x_s
    char_fit_fun[0] = 'LINFIT'
; o print result
    print, '------------------------ C_S LINFIT ------------------------'
    print, FORMAT='(A8,F10.6,A3,F5.2,A1)','coef[0]=',coeff[0],'+/-',$
           ABS(error[0]/coeff[0])*100.0,'%'
    print, FORMAT='(A8,F10.6,A3,F5.2,A1)','coef[1]=',coeff[1],'+/-',$
           ABS(error[1]/coeff[1])*100.0,'%'
    print, FORMAT='(A4,E10.3,A1,E10.3,A2,E10.3)','Cs =',CS[0,0],'+',$
           SIGMA_CS[0,0],'/-',SIGMA_CS[0,0]
    print, FORMAT='(A4,F13.6,A3,F13.6)','xs =',CS[1,0],'+/-',SIGMA_CS[2,0]
    print, 'Chi2=',CHI2_CS[0],', N=',vecsize,', prob. of fit =',probability
    print, 'covariance=',COVAB,', correlation=',RAB
    PRINT, ' no. of data points:',vecsize
    print, '------------------------ C_S LINFIT ------------------------'
ENDIF
;
fitway = 2
IF (fitway EQ 2) THEN BEGIN  ;  ------- linear regression -----------------
;   The REGRESS function performs a multiple linear regression fit and 
;   returns an Nterm-element column vector of coefficients.
;   REGRESS fits the function:
;      y_i = ln(a_i) = ln(C_s)  +  (x_s * ln(Theta_i))  +  (x_f * ln(f_i))
; o set variables
    XF      = dblarr(vecsize)
    XF[*]   = h2oh2o[0:vecsize-1,1]       ; x = ln(Theta) [1]
;
    YF      = dblarr(vecsize)
    YF[*]   = h2oh2o[0:vecsize-1,0]       ; y = ln(abs) [1]
;
    Weights = dblarr(vecsize)             
    FOR i = 0,vecsize-1 DO BEGIN
        Weights[i] = 1.00e0
    ENDFOR
;
    XSPWR98 = 4.5 ; hypothesis test with PWR98
;
; o perform the fit using multiple linear regression
    mean2sep = 1
    confa    = 1
    confy    = dblarr(N_ELEMENTS(YF))
    hypoa    = XSPWR98
    CSLR = MYLINREGRESS(XF, YF,           $
;                        weights=weights,   $ ; weight of y
                       yfit=yfit,         $ ; yfit_i = a*x_i + b
                       mean2sep=mean2sep, $ ; mean quad. separation
                       confa=confa,       $ ; 95% confidence intervall of x_s
                       confy=confy,       $ ; 95% conf.interval for each y_i
                       hypoa=hypoa)       ; x_s agreement with PWR98
;
    CS[0,1]       = EXP(CSLR[1])       ; fit values of C_s
    CS[1,1]       = CSLR[0]            ; fit values of x_s
    kk = 0
    for k = 0,N_ELEMENTS(XF)-1 do begin
        if (ABS(XF[k]) LT ABS(XF[kk])) then kk = k 
    endfor
    print,'CS  XF=',XF[kk],' confy=',confy[kk]
    if ( ABS(XF[kk]) LT 0.005 ) then begin  
        CSOp = EXP( (CSLR[1] + CSLR[0]*XF[kk]) + confy[kk] )
        CSOm = EXP( (CSLR[1] + CSLR[0]*XF[kk]) - confy[kk] )
        SIGMA_CS[0,1] = 1.000e2 * confy[kk] ; sigma_plus of C_s  in %
        SIGMA_CS[1,1] = 1.000e2 * confy[kk] ; sigma_minus of C_s in %
;        SIGMA_CS[0,1] = CS[0,1] * (EXP(confy[kk]) - 1.000e0) ; sigma_plus of C_s
;        SIGMA_CS[1,1] = CS[0,1] * (1.000e0 - EXP(confy[kk])) ; sigma_minus of C_s
    endif else begin
        SIGMA_CS[0,1] = -999.9            ; sigma_plus of C_s
        SIGMA_CS[1,1] = -999.9            ; sigma_minus of C_s
    endelse
    SIGMA_CS[2,1] = confa              ; error in x_s
    mean2sep      = mean2sep           ; sigma_plus of C_s
    char_fit_fun[0] = 'REGRESS'
; o print result
    PRINT, ' ----------------------- C_S REGRESS -----------------------'
    PRINT, 'C_s  = ',CS[0,1],' +',SIGMA_CS[0,1],' ',SIGMA_CS[1,1]
    PRINT, 'x_s  = ',CS[1,1],' +/-',SIGMA_CS[2,1]
    PRINT, 'mean2sep = ',mean2sep
    PRINT, 'hypoa_PWR98 (0=yes/1=no)=',hypoa
    PRINT, ' no. of data points:',N_ELEMENTS(YF)
    PRINT, ' ----------------------- C_S REGRESS -----------------------'
ENDIF
;
;fitway = 3
IF (fitway EQ 3) THEN BEGIN  ;  ------- simul. fit for f and Theta ---------
;   The REGRESS function performs a multiple linear regression fit and 
;   returns an Nterm-element column vector of coefficients.
;   REGRESS fits the function:
;      y_i = ln(a_i) = ln(C_s)  +  (x_s * ln(Theta_i))  +  (x_f * ln(f_i))
; o set variables
    XF       = dblarr(2,vecsize)
    XF[0,*]  = h2oh2o[0:vecsize-1,1]                ; ln(Theta) [1]
    XF[1,*]  = ALOG(h2oh2o[0:vecsize-1,2])          ; ln(f/GHz) [1]
    XF2      = REPLICATE(1.000e0, 3, vecsize)
    XF2[1,*] = h2oh2o[0:vecsize-1,1]                ; ln(Theta) [1]
    XF2[2,*] = ALOG(h2oh2o[0:vecsize-1,2])          ; ln(f/GHz) [1]
;
    YF      = dblarr(vecsize)
    YF[*]   = h2oh2o[0:vecsize-1,0]                 ; y = ln(a) [1]
;
    Weights = dblarr(vecsize)             
    FOR i = 0,vecsize-1 DO BEGIN
        Weights[i] = (1.000 / h2oh2o[i,3])^2 ; weight of each measurement
;        Weights[i] = 1.00e0
    ENDFOR
;
; o perform the fit using multiple linear regression
;    CSTf = MYREGRESS(XF, YF, Weights, yfit=yfit, Const=lnCSOTf, SIGMA=SigmaS, $
;                   FTEST=FtestS, R=RS, RMUL=RmulS, CHISQ=ChisqS)
;    CSTf = REGRESS2(XF2, YF, Weights, Yfit, SigmaS, $
;                    FtestS, RS, RmulS, ChisqS, /RELATIVE_WEIGHT)
    CS[0,2]       = EXP(CSTf[0])      ; fit values of C_s
    CS[1,2]       = CSTf[1]           ; fit values of x_s
    CS[2,2]       = CSTf[2]           ; fit values of x_s
    CHI2_CS[2]    = ChisqS            ; chi^2 of the fit
    SIGMA_CS[0,2] = CS[0,2]*(EXP( SigmaS[0])-1.0e0) ; sigma_plus of C_s
    SIGMA_CS[1,2] = CS[0,2]*(EXP(-SigmaS[0])-1.0e0) ; sigma_minus of C_s
    SIGMA_CS[2,2] = SigmaS[1]         ; error in x_s
    SIGMA_CS[3,2] = SigmaS[2]         ; error in x_f
    char_fit_fun[0] = 'REGRESS'
; o print result
    PRINT, ' ----------------------- C_S REGRESS -----------------------'
    print, 'CSTf[0]  =',CSTf[0],' (+/-',((SigmaS[0]/CSTf[0])*100.0),'%)'
    print, 'CSTf[1]  =',CSTf[1],' (+/-',((SigmaS[1]/CSTf[1])*100.0),'%)'
    print, 'CSTf[2]  =',CSTf[2],' (+/-',((SigmaS[2]/CSTf[2])*100.0),'%)'
    PRINT, ' simult. (f,Theta) fit C_s=', CS[0,2],$
           '+',SIGMA_CS[0,2],' ',SIGMA_CS[1,2]
    PRINT, ' simult. (f,Theta) fit (Theta exp.) x_s=', CS[1,2],'+/-',SIGMA_CS[2,2]
    PRINT, ' simult. (f,Theta) fit (f exp.)     y_s=', CS[2,2],'+/-',SIGMA_CS[3,2]
    PRINT, ' F for test of fit             : FtestS=',FtestS
    PRINT, ' linear correl. coeff.         : RS    =',RS
    PRINT, ' multiple linear correl. coeff.: RmulS =',RmulS
    PRINT, ' reduced, weighted chi-squared : ChisqS=',ChisqS
    PRINT, ' ----------------------- C_S REGRESS -----------------------'
ENDIF
;
;
; ==========================={MAKE A PLOT OF C_S}===========================
;
;
; --- use aii_plot_file for writing into plot output file ------------------
if (STRLEN(plotfilename) LT 3)  then plotfilename='CS_fit_plot'
if ( (plotfileformat LT 1) OR (plotfileformat GT 5) ) then plotfileformat=4
aii_plot_file, action='begin', fname=plotfilename, fformat=plotfileformat

; --- set frame of the plot in lin-lin scalo of log-log variables ----------
colors  = intarr(30) ; color array
lstyle  = intarr(30) ; line style array
pstyle  = intarr(30) ; symbol style array
aii_csf_plot_ps_definitions, colors=colors, $
  lstyle=lstyle,   $
  pstyle=pstyle,   $
  plotpos=plotpos, $
  thick=thick
;
; --- set frame of the plot ------------------------------------------------
xlnmin = MIN(h2oh2o[0:vecsize-1,1])              ; min{-ln(Theta)}_i
xlnmax = MAX(h2oh2o[0:vecsize-1,1])              ; max{-ln(Theta)}_i
ymind  = MIN(h2oh2o[0:vecsize-1,0])
ymaxd  = MAX(h2oh2o[0:vecsize-1,0])
yminm  = ALOG(CS[0,0]) + 0.20*CS[1,0] - ABS(MAX(confy[*]))
ymaxm  = ALOG(CS[0,0]) - 0.02*CS[1,0] + ABS(MAX(confy[*]))
ymin   = MIN( [ymind, yminm] )
ymax   = MAX( [ymaxd, ymaxm] )
;
plot, h2oh2o[0:vecsize-1,1],                 $ ; ln(Theta)
  h2oh2o[0:vecsize-1,0],                     $ ; ln(a)
  /NORMAL,                                   $
  xrange=[xlnmin, xlnmax],                   $
  title=TeXtoIDL('fit of H_2O-H_2O data ('+titletext+')', font=0), $
  xcharsize=1.5,                             $
  xtitle=TeXtoIDL('-ln(\Theta)     [1]', font=0),             $
  XTICK_GET = xticks,                        $
  yrange=[ymin, ymax],                       $
  ycharsize=1.5,                             $
  ytitle=TeXtoIDL('ln(C_s / [dB/km/GHz^2/hPa^2])   [1]', font=0),  $
  yTICK_GET = yticks,                        $
  color=colors[0],                           $
  psym=aii_plotsymbols(0),                   $
  symsize=1.25,                              $
  xstyle=2,                                  $
  ystyle=1,                                  $
  /nodata
;
; --- T [K] (upper x-axis) --------------------------------------------------
FOR i = 1,N_ELEMENTS(xticks)-2 DO BEGIN
    xyouts, (xticks[i]),         $ ; x coordinate
            (yticks[0]),   $ ; y coordinate
            ALIGNMENT= 0.5,      $  
            CHARSIZE= 2.0,      $   
            CHARTHICK=1.0,       $
            string((Tref*EXP(xticks[i])),  FORMAT='(I3)')+'K', $
            /DATA
ENDFOR
;
; --- plot data point for point separately ---------------------------------
fi_sym  = intarr(vecsize,2)
f_store = dblarr(vecsize)
deltax = 0.03
deltay = 0.075
k  = 0
ks = 0
FOR i = 0,vecsize-1 DO BEGIN
;   --- plot the associated legend if not done already -----
    k = find_entry_in_vec(h2oh2o[i,2], f_store[0:ks])
    IF ( k EQ -1) THEN BEGIN
        ks = ks + 1
        f_store[ks] = h2oh2o[i,2]
        fi_sym[ks,0] = pstyle[ks] ; symbols sorted by freq.
        fi_sym[ks,1] = colors[ks] ; symbol color sorted by freq.
        plots,  (plotpos[2]+deltax),                $  ; x coordinate
                (plotpos[3]-(ks*deltay)),           $  ; y coordinate
                COLOR=colors[fi_sym[ks,1]],         $
                PSYM=aii_plotsymbols(fi_sym[ks,0]), $  ; plot symbol
                SYMSIZE=1.50,                       $
                /NORMAL
        xyouts, (plotpos[2]+deltax+deltax),         $  ; x coordinate
                (plotpos[3]-0.01-(ks*deltay)),      $  ; y coordinate
                CHARSIZE= 2.0,                      $   
                CHARTHICK=1.5,                      $
                COLOR=colors[fi_sym[ks,1]],         $
                string(h2oh2o[i,2], FORMAT='(F6.2)')+' GHz', $
                /NORMAL
    ENDIF 
;   --- plot the data in the plot --------------------------
    k = find_entry_in_vec(h2oh2o[i,2], f_store[0:ks])
    plots, h2oh2o[i,1],                       $ ; ln(Theta)
           h2oh2o[i,0],                       $ ; ln(a)
           COLOR=colors[fi_sym[k,1]],         $
           PSYM=aii_plotsymbols(fi_sym[k,0]), $
           SYMSIZE=1.25
ENDFOR
;
; --- plot first fit result ------------------------------------------------ 
index = SORT(h2oh2o[*,1])
FOR j = 1,N_ELEMENTS(index)-1 DO BEGIN
    i1 = index[j-1]
    i2 = index[j]
    x1 = h2oh2o[i1,1]                         ; ln(Theta)
    x2 = h2oh2o[i2,1]                         ; ln(Theta)
    y1 = ALOG(CS[0,0]) + CS[1,0]*h2oh2o[i1,1] ; ln(Cs) + x_s*ln(Theta)
    y2 = ALOG(CS[0,0]) + CS[1,0]*h2oh2o[i2,1] ; ln(Cs) + x_s*ln(Theta)
    plots,                $
      [x1, x2],           $
      [y1, y2],           $
      color=colors[0],    $
      linestyle=0,        $
      thick=thick
ENDFOR
;
; --- plot second fit result ----------------------------------------------- 
IF (ABS(CS[1,1]) GT 0.0) THEN BEGIN
;; plot linear regression plot
    FOR j = 1,N_ELEMENTS(index)-1 DO BEGIN
        i1 = index[j-1]
        i2 = index[j]
        x1 = h2oh2o[i1,1]                           ; ln(Theta)
        x2 = h2oh2o[i2,1]                           ; ln(Theta)
        y1 = ALOG(CS[0,1]) + CS[1,1]*h2oh2o[i1,1]   ; ln(Cs) + x_s*ln(Theta)
        y2 = ALOG(CS[0,1]) + CS[1,1]*h2oh2o[i2,1]   ; ln(Cs) + x_s*ln(Theta)        
        dy1  = ALOG(CS[0,1]) + CS[1,1]*h2oh2o[i1,1] - confy[i1]
        dy2  = ALOG(CS[0,1]) + CS[1,1]*h2oh2o[i2,1] - confy[i2]
        ddy1 = ALOG(CS[0,1]) + CS[1,1]*h2oh2o[i1,1] + confy[i1]
        ddy2 = ALOG(CS[0,1]) + CS[1,1]*h2oh2o[i2,1] + confy[i2]
        plots, $
          [x1, x2],           $ ; ln(Theta)
          [y1, y2],           $ ; ln(Cs) + x_s*ln(Theta)
          color=colors[6],    $
          linestyle=1,        $
          thick=thick
        if ( ( (dy1 GT ymin) AND (dy1 LT ymax) ) AND $
             ( (dy2 GT ymin) AND (dy2 LT ymax) ) ) then begin
            plots, $
              [x1, x2],           $ ; ln(Theta)
              [dy1, dy2],         $ ; ln(Cs) + x_s*ln(Theta)
              color=colors[6],    $
              linestyle=1,        $
              thick=thick
        endif
        if ( ( (ddy1 GT ymin) AND (ddy1 LT ymax) ) AND $
             ( (ddy2 GT ymin) AND (ddy2 LT ymax) ) ) then begin
            plots, $
              [x1, x2],           $ ; ln(Theta)
              [ddy1, ddy2],       $ ; ln(Cs) + x_s*ln(Theta)
              color=colors[6],    $
              linestyle=1,        $
              thick=thick
        endif
    ENDFOR
ENDIF
;
; --- write date and user info ---------------------------------------------
datum = get_userdate_info(1)
xyouts, plotpos[0], plotpos[1]-0.10, datum, CHARSIZE=0.75, CHARTHICK=1.0, /NORMAL
;
; close plot output file
aii_plot_file, action='end', show='no', print='no', $
               outdir='./' 
;
;goto, omit_print
;
; =========================={ PRINT THE STRUCTURE }=========================
;
;
openw, xunit, artsjobpath+controlfilenamecore+'.param', ERROR=err, /get_lun
IF (err NE 0) THEN BEGIN
    print,' do_PWR98_cont_fit> ERROR: output file >>'+artsjobpath+controlfilenamecore+$
          '.param'+'<< can not be created!'
    print,' do_PWR98_cont_fit> error code for output file:',err
    print,' do_PWR98_cont_fit> error message:', !ERR_STRING
    FREE_LUN, xunit
ENDIF
printf, xunit,'# =========================================================================='
printf, xunit,'# Data to fit the water vapor continuum parameters Cs, xs, Cf, xf'
printf, xunit,'# row entry information:'
printf, xunit,'# 0: buffer gas species '
printf, xunit,'# 1: buffer gas pressure              [Pa]'
printf, xunit,'# 2: H2O species'
printf, xunit,'# 3: H2O pressure                     [Pa]'
printf, xunit,'# 4: frequency of measurement         [GHz]'
printf, xunit,'# 5: termperature of measurement      [K]'
printf, xunit,'# 6: measured total absorption        [1/m]'
printf, xunit,'# 7: line absorption model name'
printf, xunit,'# 8: calculated line absorption       [1/m]'
printf, xunit,'# 9: calculated continuum absorption  [1/m]'
printf, xunit,'# =========================================================================='
FOR i = 0, N_ELEMENTS(data[*].Abs)-1 DO BEGIN
    printf, xunit,  $
           FORMAT='(I3,A1,A3,A1,F8.1,A1,A3,A1,F6.1,A1,F7.2,A1,F5.1,A1,E10.3,A1,A20,A1,E10.3,A1,E10.3)', $
           i,' ', $
           data[i].BuffName, ' ',$
           data[i].pd, ' ', $
           data[i].WVName, ' ', $
           data[i].pwv, ' ', $
           data[i].f*1.0e-9, ' ', $
           data[i].T, ' ', $
           data[i].abs, ' ', $
           data[i].LineMod, ' ', $
           data[i].labs, ' ', $
           data[i].cabs
ENDFOR
FREE_LUN, xunit
;
; ---- END OF FUNCTION -------------------------------------
omit_print:
read_lineabs_data_ende:
; set saved settings back
!P = P_ini
IF (CS[0,1] LE 0.0) THEN BEGIN
; [C_s, sigma_C_s, x_s, sigma_x_s, chi2, cov(a,b), correlation(a,b)]
    resvec = dblarr(7)          
    resvec[0] = CS[0,0]
    resvec[1] = 0.5 * ( SIGMA_CS[0,0] + SIGMA_CS[1,0] )
    resvec[2] = CS[1,0]
    resvec[3] = SIGMA_CS[2,0]
    resvec[4] = CHI2_CS[0]  
    resvec[5] = COVAB
    resvec[6] = RAB
ENDIF ELSE BEGIN
    resvec = dblarr(12)
    resvec[0]  = CS[0,0]
    resvec[1]  = 0.5 * ( SIGMA_CS[0,0] + SIGMA_CS[1,0] )
    resvec[2]  = CS[1,0]
    resvec[3]  = SIGMA_CS[2,0]
    resvec[4]  = CHI2_CS[0]  
    resvec[5]  = COVAB
    resvec[6]  = RAB
    resvec[7]  = CS[0,1]       ; linear regression fit
    resvec[8]  = CS[1,1]       ; linear regression fit
    resvec[9]  = SIGMA_CS[0,1] ; linear regression fit
    resvec[10] = SIGMA_CS[1,1] ; linear regression fit
    resvec[11] = SIGMA_CS[2,1] ; linear regression fit
ENDELSE

RETURN,  resvec
END
;
;
; **************************************************************************
; NAME:     do_PWR98_h2on2_cont_fit
;
; PURPOSE:  Perform a H2O continuum parameter fit according to the water
;           vapor continuum parameterization of P. W. Rosenkranz,
;           Radio Science, 33(4),  919, 1998, Radio Science, 34(4), 1025, 1999
;           The only data considered is H2O-N2 mixture.
;
; INPUTS:   filename      idl structure with the all data
;
; OUTPUT:   flag          0: ok  / else: false
;           file          the file is named after the core file name
;                         of this job with the extension '.param'
;
; EXTERNAL: TeXtoIDL package (http://physweb.mnstate.edu/mcraig/TextoIDL/v2-0/default.htm)
;           CreateArtsControlFile.pro (arts/aii)
;           aii_checks.pro            (arts/aii)
;           aii_plot_file.pro         (arts/aii)
;           aii_plotsymbols.pro       (arts/aii)
;
; HISTORY:  2002-01-16    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION do_PWR98_h2on2_cont_fit, data,                             $
                                  CS, sigma_Cs, xs, sigma_xs,       $
                                  artsjobpath, controlfilenamecore, $
                                  linetotmaxratio,                  $
                                  FFITMIN, FFITMAX,                 $
                                  plotfilename, plotfileformat,     $
                                  titletext
;
;
; ---- CHECK OF CS ---------------------------------------------------------
IF (CS LE 0.0000e0) THEN BEGIN 
    print, ' do_PWR98_h2on2_cont_fit> ERROR detected!'
    print, ' do_PWR98_h2on2_cont_fit> wrong input parameter: C_s=',CS
ENDIF
IF (xs LE -1.0000e1) THEN BEGIN 
    print, ' do_PWR98_h2on2_cont_fit> ERROR detected!'
    print, ' do_PWR98_h2on2_cont_fit> wrong input parameter: x_s=',xs
ENDIF
IF (sigma_cs LE 0.0000e0) THEN BEGIN 
    print, ' do_PWR98_h2on2_cont_fit> ERROR detected!'
    print, ' do_PWR98_h2on2_cont_fit> wrong input parameter: sigma_C_s=',sigma_cs
ENDIF
IF (sigma_xs LE 0.0000e0) THEN BEGIN 
    print, ' do_PWR98_h2on2_cont_fit> ERROR detected!'
    print, ' do_PWR98_h2on2_cont_fit> wrong input parameter: sigma_x_s=',sigma_xs
ENDIF
print,' do_PWR98_h2on2_cont_fit> 0 CS=',CS,'(',sigma_Cs,') xs=', xs,'(',sigma_xs,')'
;CS       = 9.111e-08
;sigma_Cs = 6.751E-09
xs       = -xs
;sigma_xs = 0.190
print,' do_PWR98_h2on2_cont_fit> 1 CS=',CS,'(',sigma_Cs,') xs=', xs,'(',sigma_xs,')'
;
; ---- OUTPUT VARIABLE -----------------------------------------------------
flag = 1
;
; ---- SAVE SETTINGS -------------------------------------------------------
P_ini = !P
;
; ---- GENERAL CONSTANTS ---------------------------------------------------
COMMON UNITCONVERSION, TENLOG10_EULER, dBkm2Npm, Npm2dBkm, Hz2GHz, Pa2hPa
;
Tref                 = 3.000E2    ; [K] for Theta calculation
Rlinetotmax          = 1.00000    ; ratio of line to total absorption
Rlinetotmin          = 0.01       ; ratio of line to total absorption
datanlim             = 5          ; minimum of data points for fit
sigma_ratio_abs_tot  = 1.000e-1   ; Gaussian error propagation: sigma_abstot = 10.0% abs_tot 
sigma_ratio_abs_line = 5.000e-2   ; Gaussian error propagation: sigma_absl   =  5.0% abs_l
sigma_ratio_f        = 1.000e-5   ; Gaussian error propagation: sigma_f      =  0.001% f
sigma_ratio_pwv      = 1.000e-3   ; Gaussian error propagation: sigma_PH2O   =  0.1% P_H2O
sigma_ratio_pd       = 1.000e-4   ; Gaussian error propagation: sigma_buffer =  0.01% P_buffer
sigma_ratio_T        = 5.000e-3   ; Gaussian error propagation. sigma_T      =  0.5% T
IF (linetotmaxratio GE 1.00) THEN linetotmaxratio = Rlinetotmax 
;
; --- FIT VARIABLES H2O-H2O ------------------------------------------------
; these are the variables with the fit results for max. n_fits different
; fit procedures
n_fits       = 3
CF           = dblarr(3,n_fits) ; C_s for different fit routines
                                ; 0 : C_s 
                                ; 1 : x_s 
                                ; 2 : x_f (if performed) 
CF[*,*]      = 0.000e0
SIGMA_CF     = dblarr(4,n_fits) ; std. deviations for different fit routines
                                ; 0 : +sigma of C_s
                                ; 1 : -sigma of C_s
                                ; 2 : sigma of x_s
                                ; 3 : sigma of x_f (if freq. fit is performed)
SIGMA_CF[*,*]= 0.000e0
CHI2_CF      = dblarr(n_fits)   ; chi² of the different fit routines
CHI2_CF[*]   = 0.000e0
COVAB = 0.000
RAB   = 0.000
char_fit_fun = strarr(n_fits)   ; strings, names the for different fit routines
;
; --- SORT N2-H2O DATA -----------------------------------------------------
h2on2 = dblarr(N_ELEMENTS(data[*].abs), 4) ; 0 : y =  ln(a)
                                           ; 1 : x = -ln(Theta)
                                           ; 2 : frequency [GHz]
                                           ; 3 : sigma_y error estimation of y
no = 0
j  = 0
FOR i = 0, N_ELEMENTS(data[*].abs)-1 DO BEGIN
;    print, data[i].abs,'|',data[i].labs,'|',data[i].BuffName ,'|',data[i].WVName
;    print, data[i].f,'|',data[i].T,'|',data[i].pd,'|',data[i].pwv
;    print, '--------------------------------------------------------------------'
    IF ((data[i].BuffName EQ 'N2') AND (data[i].WVName EQ 'H2O'))  THEN BEGIN
;        print,' do_PWR98_cont_fit> f=',data[i].f,' r=',(data[i].labs / data[i].abs)
;   --- check frequency range ------------------------------------------------
        IF ( ( data[i].f GE FFITMIN ) AND $
             ( data[i].f LE FFITMAX )   ) THEN BEGIN
;   --- check line to total absorption ratio ---------------------------------
        IF ( ((data[i].labs / data[i].abs) LT linetotmaxratio) AND $
             ((data[i].labs / data[i].abs) GT Rlinetotmin) ) THEN BEGIN
;         o Theta [1]
            Theta = (Tref / data[i].T)
;         o effective "absorption" term [dB/km]
            z = Npm2dBkm * ( data[i].abs - data[i].labs )
;         o nominator [GHz²*hPa²]
            n = (Hz2GHz*data[i].f)^2    *  $ ; f²     [GHz²]
                Theta^3                 *  $ ; Theta³ [1]
                (Pa2hPa*data[i].pwv)    *  $ ; p_H2O  [hPa]
                (Pa2hPa*data[i].pd)          ; p_N2   [hPa]
            a = (z/n)                        ; => a_i = (z/n)
;         o self contribution = C_s * Theta^x-s * P_H2O/P_N2  [dB/km/GHz²/hPa²] 
            s = CS * Theta^xs * (data[i].pwv / data[i].pd)
;         o effective foreign contribution to fit
            IF (( a - s )  GT 0.0e0) THEN BEGIN
;             o "effective absorption" to fit
                h2on2[j,0] =  ALOG( a - s ) ; => y_i = ln(a_i) with a_i = (z/n)
;             o Ln(Theta) [1]
                h2on2[j,1] = -ALOG( Theta ) ; => x_i = -ln(Theta)
;             o frequency [Hz]
                h2on2[j,2] = Hz2GHz * data[i].f
;             o calculate an error estimate for  absorption [1/m/GHz²/hPa²]
;               with Gaussian error propagation
                diffa    = dblarr(8)
                siga     = dblarr(8)
                diffa[0] =  1.000e0/n                              ; d a / d abstot
                diffa[1] = -1.000e0/n                              ; d a / d absl
                diffa[2] =  a * (-2.00e0/(Hz2GHz*data[i].f))       ; d a / d f 
                diffa[3] =  a * (-1.00e0/(Pa2hPa*data[i].pwv)) - $ ; d a / d P_H2O
                            (s/(Pa2hPa*data[i].pwv))                
                diffa[4] =  a * (-1.00e0/(Pa2hPa*data[i].pd))  + $ ; d a / d P_N2
                            (s/(Pa2hPa*data[i].pd))                  
                diffa[5] =  a * (-3.00e0/Theta) - s * (xs/Theta)   ; d a / d Theta
                diffa[6] = s / CS                                  ; d a / d Cs
                diffa[7] = s * ALOG(Theta)                         ; d a / d xs
                siga[0]  = sigma_ratio_abs_tot  * (Npm2dBkm*data[i].abs)
                siga[1]  = sigma_ratio_abs_line * (Npm2dBkm*data[i].labs)
                siga[2]  = sigma_ratio_f        * (Hz2GHz*data[i].f)
                siga[3]  = sigma_ratio_pwv      * (Pa2hPa*data[i].pwv)
                siga[4]  = sigma_ratio_pd       * (Pa2hPa*data[i].pd)
                siga[5]  = sigma_ratio_T        * Theta
                siga[6]  = sigma_Cs
                siga[7]  = sigma_xs
                sig_a    = 0.000e0
                FOR k = 0, N_ELEMENTS(diffa)-1 DO BEGIN
                    sig_a = sig_a + ( (diffa[k])^2 * (siga[k])^2 ) 
;                    print, k,': diffa=',diffa[k],', siga=',siga[k],', sig_a=',sig_a
                ENDFOR
                sig_a      = SQRT( sig_a )
                sig_y_p    = ABS( ALOG((z/n)+sig_a) - h2on2[j,0] ) ; upper side
                sig_y_m    = ABS( ALOG((z/n)-sig_a) - h2on2[j,0] ) ; lower side
                h2on2[j,3] = 0.500e0*(sig_y_p+sig_y_m)             ; take the mean of both
;             o print data for control
;                print,j,': x=',h2on2[j,1],', y=',h2on2[j,0],', s_y=',h2on2[j,3]
;             o set data counter one higher
                j = j + 1
            ENDIF ELSE BEGIN
                no = no + 1
            ENDELSE
        ENDIF
        ENDIF
    ENDIF
ENDFOR
print,' do_PWR98_cont_fit> # of H2O-N2 data points omitted=',no
IF (j LT datanlim) THEN BEGIN
    print,' do_PWR98_cont_fit> too few data points (#=',j,') for H2O-N2 fit.'
    print,' do_PWR98_cont_fit> jump to the end of this function!'
    goto, read_lineabs_data_ende
ENDIF
print,' do_PWR98_cont_fit> no. of H2O-N2 data points for fit=',j
h2on22  = h2on2[0:j-1, 0:N_ELEMENTS(h2on2[0,*])-1] ; resize array
index   = SORT(h2on22[*,1])                        ; sort in increasing ln(Theta)
vecsize = N_ELEMENTS(h2on22[*,0])                  ; just for simplicity
IF (N_ELEMENTS(index) NE vecsize) THEN BEGIN
    print,' do_PWR98_cont_fit> !!! ERROR: inconsistency in vector size!'
    print,' do_PWR98_cont_fit> !!! N_ELEMENTS(index):',N_ELEMENTS(index)
    print,' do_PWR98_cont_fit> !!! vecsize          :',vecsize
    goto, read_lineabs_data_ende
ENDIF
FOR i = 0, vecsize-1 DO BEGIN
    j = index[i]
    h2on2[i,0] = h2on22[j,0] ; order the array in increasing ln(Theta)
    h2on2[i,1] = h2on22[j,1]
    h2on2[i,2] = h2on22[j,2]
    h2on2[i,3] = h2on22[j,3]
    print, FORMAT='(I2,A4,F8.4,A6,F8.4,A4,F6.3,A4,F7.3)', $
           i,': y=',h2on2[i,0],': s_y=',h2on2[i,3],', x=',h2on2[i,1],', f=',h2on2[i,2]
ENDFOR
h2on2 = h2on22[0:vecsize-1, 0:N_ELEMENTS(h2on2[0,*])-1]      ; resize array
;
;
; ==========================={FIT H2O-N2 DATA}==============================
;
;
fitway = 1
IF (fitway EQ 1) THEN BEGIN ;  ----------- LINFIT --------------------------
;   The LINFIT function fits the paired data {xi, yi} to the 
;   linear model, y = A + Bx, by minimizing the Chi-square error statistic. 
;   The result is a two-element vector containing the model parameters [A, B]. 
; o perform fit
    sdevvec = dblarr( N_ELEMENTS(h2on2[0:vecsize-1,0]) )
    for isdev = 0, N_ELEMENTS(sdevvec)-1 do $
      sdevvec[isdev] = 1.000 / ABS( h2on2[isdev,0] )
    coeff = MYLINFIT(h2on2[0:vecsize-1,1],    $
                   h2on2[0:vecsize-1,0],      $
;                   SDEV=sdevvec,              $
                   PROB=probability , $
                   SIGMA=error,       $
                   CHISQ=CHISQ2,      $
                   COVAB=COVAB,       $
                   RAB=RAB,           $
                   /DOUBLE)
    CF[0,0]       = EXP(coeff[0])                   ; fit values of C_s
    CF[1,0]       = coeff[1]                        ; fit values of x_s
    CHI2_CF[0]    = CHISQ2                          ; chi^2 of the fit
    SIGMA_CF[0,0] = CF[0,0]*(EXP( error[0])-1.0e0)  ; sigma_plus of C_s
    SIGMA_CF[1,0] = CF[0,0]*(EXP(-error[0])-1.0e0)  ; sigma_minus of C_s
    SIGMA_CF[2,0] = error[1]                        ; error in x_s
    char_fit_fun[0] = 'LINFIT'
; o print result
    print, '------------------------ C_F LINFIT ------------------------'
    print, FORMAT='(A8,F10.6,A3,F5.2,A1)','coef[0]=',coeff[0],'+/-',$
           ABS(error[0]/coeff[0])*100.0,'%'
    print, FORMAT='(A8,F10.6,A3,F5.2,A1)','coef[1]=',coeff[1],'+/-',$
           ABS(error[1]/coeff[1])*100.0,'%'
    print, FORMAT='(A4,E10.3,A3,E10.3,A3,E10.3)','Cf =',CF[0,0],' /+',$
           SIGMA_CF[0,0],' /-',SIGMA_CF[0,0]
    print, FORMAT='(A4,F13.6,A3,F13.6)','xf =',CF[1,0],'+/-',SIGMA_CF[2,0]
    print, 'Chi2=',CHI2_CF[0],', N=',vecsize,', prob. of fit =',probability
    print, 'covariance=',COVAB,', correlation=',RAB
    PRINT, ' no. of data points:',vecsize
    print, '------------------------ C_F LINFIT ------------------------'
ENDIF
;
;
fitway = 2
IF (fitway EQ 2) THEN BEGIN  ;  ------- linear regression -----------------
;   The REGRESS function performs a multiple linear regression fit and 
;   returns an Nterm-element column vector of coefficients.
;   REGRESS fits the function:
;      y_i = ln(a_i) = ln(C_s)  +  (x_s * ln(Theta_i))  +  (x_f * ln(f_i))
; o set variables
    XF      = dblarr(vecsize)
    XF[*]   = h2on2[0:vecsize-1,1]       ; x = ln(Theta) [1]
;
    YF      = dblarr(vecsize)
    YF[*]   = h2on2[0:vecsize-1,0]       ; y = ln(abs) [1]
;
    Weights = dblarr(vecsize)             
    FOR i = 0,vecsize-1 DO BEGIN
        Weights[i] = 1.00e0
    ENDFOR
;
    XSPWR98 = 0.0 ; hypothesis test with PWR98
;
; o perform the fit using multiple linear regression
    mean2sep = 1
    confa    = 1
    confy    = dblarr(N_ELEMENTS(YF))
    hypoa    = XSPWR98
    CFLR = MYLINREGRESS(XF, YF,           $
;                        weights=weights,   $ ; weight of y
                       yfit=yfit,         $ ; yfit_i = a*x_i + b
                       mean2sep=mean2sep, $ ; mean quad. separation
                       confa=confa,       $ ; 95% confidence intervall of a
                       confy=confy,       $ ; 95% conf.interval for each y_i
                       hypoa=hypoa)       ; x_s agreement with PWR98
;
    CF[0,1]       = EXP(CFLR[1])       ; fit values of C_f
    CF[1,1]       = CFLR[0]            ; fit values of x_f
    kk = 0
    for k = 0,N_ELEMENTS(XF)-1 do begin
        if (ABS(XF[k]) LT ABS(XF[kk])) then kk = k 
    endfor
    print,'CF  XF=',XF[kk],' confy=',confy[kk]
    if ( ABS(XF[kk]) LT 0.005 ) then begin  
        CFOp = EXP( (CFLR[1] + CFLR[0]*XF[kk]) + confy[kk] )
        CFOm = EXP( (CFLR[1] + CFLR[0]*XF[kk]) - confy[kk] )
        SIGMA_CF[0,1] = 1.000e2 * confy[kk] ; sigma_plus of C_f  in %
        SIGMA_CF[1,1] = 1.000e2 * confy[kk] ; sigma_minus of C_f in %
    endif else begin
        SIGMA_CF[0,1] = -999.9            ; sigma_plus of C_f
        SIGMA_CF[1,1] = -999.9            ; sigma_minus of C_f
    endelse
    SIGMA_CF[2,1] = confa              ; error in x_f
    mean2sep      = mean2sep           ; sigma_plus of C_s
    char_fit_fun[0] = 'REGRESS'
; o print result
    PRINT, ' ----------------------- C_F REGRESS -----------------------'
    PRINT, 'C_f  = ',CF[0,1],' +',SIGMA_CF[0,1],' -',SIGMA_CF[1,1]
    PRINT, 'x_f  = ',CF[1,1],' +/-',SIGMA_CF[2,1]
    PRINT, 'mean2sep = ',mean2sep
    PRINT, 'hypoa_PWR98 (0=yes/1=no)=',hypoa
    PRINT, ' no. of data points:',N_ELEMENTS(YF)
    PRINT, ' ----------------------- C_F REGRESS -----------------------'
ENDIF
;
;fitway = 3
IF (fitway EQ 3) THEN BEGIN  ;  ------- simul. fit for f and Theta ---------
;   The REGRESS function performs a multiple linear regression fit and 
;   returns an Nterm-element column vector of coefficients.
;   REGRESS fits the function:
;      y_i = ln(a_i) = ln(C_s)  +  (x_s * ln(Theta_i))  +  (x_f * ln(f_i))
; o set variables
    XF       = dblarr(2,vecsize)
    XF[0,*]  = h2on2[0:vecsize-1,1]                ; ln(Theta) [1]
    XF[1,*]  = ALOG(h2on2[0:vecsize-1,2])          ; ln(f/GHz) [1]
    XF2      = REPLICATE(1.000e0, 3, vecsize)
    XF2[1,*] = h2on2[0:vecsize-1,1]                ; ln(Theta) [1]
    XF2[2,*] = ALOG(h2on2[0:vecsize-1,2])          ; ln(f/GHz) [1]
;
    YF      = dblarr(vecsize)
    YF[*]   = h2on2[0:vecsize-1,0]                 ; y = ln(a) [1]
;
    Weights = dblarr(vecsize)             
    FOR i = 0,vecsize-1 DO BEGIN
        Weights[i] = (1.000 / h2on2[i,3])^2 ; weight of each measurement
;        Weights[i] = 1.00e0
    ENDFOR
;
; o perform the fit using multiple linear regression
    CFTf = MYREGRESS(XF, YF, Weights, yfit=yfit, Const=lnCSOTf, Sigma=SigmaS, $
                     FTEST=FtestS, R=RS, RMUL=RmulS, CHISQ=ChisqS)
;    CFTf = REGRESS2(XF2, YF, Weights, Yfit, SigmaS, $
;                    FtestS, RS, RmulS, ChisqS, RELATIVE_WEIGHT=0)
    CF[0,2]       = EXP(CFTf[0])      ; fit values of C_s
    CF[1,2]       = CFTf[1]           ; fit values of x_s
    CF[2,2]       = CFTf[2]           ; fit values of x_s
    CHI2_CF[2]    = ChisqS            ; chi^2 of the fit
    SIGMA_CF[0,2] = CF[0,2]*(EXP( SigmaS[0])-1.0e0) ; sigma_plus of C_s
    SIGMA_CF[1,2] = CF[0,2]*(EXP(-SigmaS[0])-1.0e0) ; sigma_minus of C_s
    SIGMA_CF[2,2] = SigmaS[1]         ; error in x_s
    SIGMA_CF[3,2] = SigmaS[2]         ; error in x_f
    char_fit_fun[0] = 'REGRESS'
    COVAB = 0.00
    RAB   = 0.00
; o print result
    PRINT, ' ----------------------- C_F REGRESS -----------------------'
    print, 'CFTf[0]  =',CFTf[0],' (+/-',((SigmaS[0]/CFTf[0])*100.0),'%)'
    print, 'CFTf[1]  =',CFTf[1],' (+/-',((SigmaS[1]/CFTf[1])*100.0),'%)'
    print, 'CFTf[2]  =',CFTf[2],' (+/-',((SigmaS[2]/CFTf[2])*100.0),'%)'
    PRINT, ' simult. (f,Theta) fit C_f=', CF[0,2],$
           '+',SIGMA_CF[0,2],' ',SIGMA_CF[1,2]
    PRINT, ' simult. (f,Theta) fit (Theta exp.) x_f=', CF[1,2],'+/-',SIGMA_CF[2,2]
    PRINT, ' simult. (f,Theta) fit (f exp.)     y_f=', CF[2,2],'+/-',SIGMA_CF[3,2]
    PRINT, ' F for test of fit             : Ftest =',FtestS
    PRINT, ' linear correl. coeff.         : R     =',RS
    PRINT, ' multiple linear correl. coeff.: Rmul  =',RmulS
    PRINT, ' reduced, weighted chi-squared : Chisq =',ChisqS
    PRINT, ' ----------------------- C_F REGRESS -----------------------'
ENDIF
;
;
;
; ==========================={MAKE A PLOT OF C_F}===========================
;
;
; --- set plot environment -------------------------------------------------
if (STRLEN(plotfilename) LT 3)  then plotfilename='CF_fit_plot'
if ( (plotfileformat LT 1) OR (plotfileformat GT 5) ) then plotfileformat=4
aii_plot_file, action='begin', fname=plotfilename, fformat=plotfileformat

; --- set frame of the plot in lin-lin scalo of log-log variables ----------
colors  = intarr(30) ; color array
lstyle  = intarr(30) ; line style array
pstyle  = intarr(30) ; symbol style array
aii_csf_plot_ps_definitions, colors=colors, $
  lstyle=lstyle,   $
  pstyle=pstyle,   $
  plotpos=plotpos, $
  thick=thick
;;
; --- set frame of the plot ------------------------------------------------
xlnmin = MIN(h2on2[0:vecsize-1,1])              ; min{-ln(Theta)}_i
xlnmax = MAX(h2on2[0:vecsize-1,1])              ; max{-ln(Theta)}_i
ymind  = MIN(h2on2[0:vecsize-1,0])
ymaxd  = MAX(h2on2[0:vecsize-1,0])
yminm  = ALOG(CF[0,0]) + 0.20*CF[1,0] - ABS(MAX(confy[*]))
ymaxm  = ALOG(CF[0,0]) - 0.02*CF[1,0] + ABS(MAX(confy[*]))
ymin   = MIN( [ymind, yminm] )
ymax   = MAX( [ymaxd, ymaxm] )
;
plot, h2on2[0:vecsize-1,1],                  $  ; x=-ln(Theta)
  h2on2[0:vecsize-1,0],                      $  ; y= ln(a)
  /NORMAL,                                   $
  xrange=[xlnmin, xlnmax],                   $
  title=TeXtoIDL('fit of H_2O-N_2 data ('+titletext+')', font=0), $
  xcharsize=1.5,                             $
  xtitle=TeXtoIDL('-ln(\Theta)     [1]', font=0),                    $
  XTICK_GET = xticks,                        $
  yrange=[ymin, ymax],                       $
  ycharsize=1.5,         $
  ytitle=TeXtoIDL('ln(C_f / [dB/km/GHz^2/hPa^2])   [1]', font=0), $
  yTICK_GET = yticks,                        $
  color=colors[0],                           $
  psym=aii_plotsymbols(0),                   $
  symsize=1.25,                              $
  xstyle=2,                                  $
  ystyle=2,                                  $
  /nodata
;
; --- T [K] (upper x-axis) -------------------------------------------------
FOR i = 1,N_ELEMENTS(xticks)-2 DO BEGIN
    xyouts, (xticks[i]),         $ ; x coordinate
            (yticks[0]),         $ ; y coordinate
            ALIGNMENT= 0.5,      $  
            CHARSIZE= 2.0,       $   
            CHARTHICK=1.0,       $
            string((Tref*EXP(xticks[i])),  FORMAT='(I3)')+'K', $
            /DATA
ENDFOR
;
; --- plot data point for point separately ---------------------------------
fi_sym  = intarr(vecsize,2)
f_store = dblarr(vecsize)
deltax = 0.03
deltay = 0.075
k  = 0
ks = 0
FOR i = 0,vecsize-1 DO BEGIN
;   --- plot the associated legend if not done already -----
    k = find_entry_in_vec(h2on2[i,2], f_store[0:ks])
    IF ( k EQ -1) THEN BEGIN
        ks = ks + 1
        f_store[ks] = h2on2[i,2]
        fi_sym[ks,0] = pstyle[ks] ; symbols sorted by freq.
        fi_sym[ks,1] = colors[ks] ; symbol color sorted by freq.
        plots,  (plotpos[2]+deltax),                $  ; x coordinate
                (plotpos[3]-(ks*deltay)),           $  ; y coordinate
                COLOR=colors[fi_sym[ks,1]],         $
                PSYM=aii_plotsymbols(fi_sym[ks,0]), $  ; plot symbol
                SYMSIZE=1.50,                       $
                /NORMAL
        xyouts, (plotpos[2]+deltax+deltax),         $  ; x coordinate
                (plotpos[3]-0.01-(ks*deltay)),      $  ; y coordinate
                CHARSIZE= 2.0,                      $   
                CHARTHICK=1.5,                      $
                COLOR=colors[fi_sym[ks,1]],         $
                string(h2on2[i,2], FORMAT='(F6.2)')+' GHz', $
                /NORMAL
    ENDIF 
;   --- plot the data in the plot --------------------------
    k = find_entry_in_vec(h2on2[i,2], f_store[0:ks])
    plots, h2on2[i,1],                        $ ; ln(Theta)
           h2on2[i,0],                        $ ; ln(a)
           COLOR=colors[fi_sym[k,1]],         $
           PSYM=aii_plotsymbols(fi_sym[k,0]), $
           SYMSIZE=1.25
ENDFOR
;
; --- plot first fit result ------------------------------------------------ 
index = SORT(h2on2[*,1])
FOR j = 1,N_ELEMENTS(index)-1 DO BEGIN
    i1 = index[j-1]
    i2 = index[j]
    x1 = h2on2[i1,1]                         ; ln(Theta)
    x2 = h2on2[i2,1]                         ; ln(Theta)
    y1 = ALOG(CF[0,0]) + CF[1,0]*h2on2[i1,1] ; ln(Cs) + x_s*ln(Theta)
    y2 = ALOG(CF[0,0]) + CF[1,0]*h2on2[i2,1] ; ln(Cs) + x_s*ln(Theta)
    plots,                $
      [x1, x2],           $
      [y1, y2],           $  
      color=colors[0],    $
      linestyle=0,        $
      thick=thick
ENDFOR
;
; --- plot second fit result ----------------------------------------------- 
IF (ABS(CF[1,1]) GT 0.0) THEN BEGIN
;; plot linear regression plot
    FOR j = 1,N_ELEMENTS(index)-1 DO BEGIN
        i1 = index[j-1]
        i2 = index[j]
        x1 = h2on2[i1,1]                           ; ln(Theta)
        x2 = h2on2[i2,1]                           ; ln(Theta)
        y1 = ALOG(CF[0,1]) + CF[1,1]*h2on2[i1,1]   ; ln(Cs) + x_s*ln(Theta)
        y2 = ALOG(CF[0,1]) + CF[1,1]*h2on2[i2,1]   ; ln(Cs) + x_s*ln(Theta)        
        dy1  = ALOG(CF[0,1]) + CF[1,1]*h2on2[i1,1] - confy[i1]
        dy2  = ALOG(CF[0,1]) + CF[1,1]*h2on2[i2,1] - confy[i2]
        ddy1 = ALOG(CF[0,1]) + CF[1,1]*h2on2[i1,1] + confy[i1]
        ddy2 = ALOG(CF[0,1]) + CF[1,1]*h2on2[i2,1] + confy[i2]
        plots, $
          [x1, x2],           $ ; ln(Theta)
          [y1, y2],           $ ; ln(Cs) + x_s*ln(Theta)
          color=colors[6],    $
          linestyle=1,        $
          thick=thick
        if ( ( (dy1 GT ymin) AND (dy1 LT ymax) ) AND $
             ( (dy2 GT ymin) AND (dy2 LT ymax) ) ) then begin
            plots, $
              [x1, x2],           $ ; ln(Theta)
              [dy1, dy2],         $ ; ln(Cs) + x_s*ln(Theta)
              color=colors[6],    $
              linestyle=1,        $
              thick=thick
        endif
        if ( ( (ddy1 GT ymin) AND (ddy1 LT ymax) ) AND $
             ( (ddy2 GT ymin) AND (ddy2 LT ymax) ) ) then begin
            plots, $
              [x1, x2],           $ ; ln(Theta)
              [ddy1, ddy2],       $ ; ln(Cs) + x_s*ln(Theta)
              color=colors[6],    $
              linestyle=1,        $
              thick=thick
        endif
    ENDFOR
ENDIF
;
; --- write date and user info ---------------------------------------------
datum = get_userdate_info(1)
xyouts, plotpos[0], plotpos[1]-0.10, datum, CHARSIZE=0.75, CHARTHICK=1.0, /NORMAL
;
; close plot output file
aii_plot_file, action='end', show='no', print='no', $
               outdir='./' 
;
; ---- END OF FUNCTION -------------------------------------
read_lineabs_data_ende:
; set saved settings back
!P = P_ini
IF (CF[0,1] LE 0.00) THEN BEGIN
; result vector [C_f, sigma_C_f, x_f, sigma_x_f, chi2, cov(a,b), correlation(a,b)]
    resvec = dblarr(7) 
    resvec[0]  = CF[0,0]
    resvec[1]  = 0.5 * ( SIGMA_CF[0,0] + SIGMA_CF[1,0] ) ; take mean as sigma
    resvec[2]  = CF[1,0]
    resvec[3]  = SIGMA_CF[2,0]
    resvec[4]  = CHI2_CF[0]  
    resvec[5]  = COVAB
    resvec[6]  = RAB
ENDIF ELSE BEGIN
    resvec = dblarr(12)
    resvec[0]  = CF[0,0]
    resvec[1]  = 0.5 * ( SIGMA_CF[0,0] + SIGMA_CF[1,0] ) ; take mean as sigma
    resvec[2]  = CF[1,0]
    resvec[3]  = SIGMA_CF[2,0]
    resvec[4]  = CHI2_CF[0]  
    resvec[5]  = COVAB
    resvec[6]  = RAB
    resvec[7]  = CF[0,1]       ; linear regression fit
    resvec[8]  = CF[1,1]       ; linear regression fit
    resvec[9]  = SIGMA_CF[0,1] ; linear regression fit
    resvec[10] = SIGMA_CF[1,1] ; linear regression fit
    resvec[11] = SIGMA_CF[2,1] ; linear regression fit
ENDELSE
RETURN, resvec
END
;
;
;
;############################# MAIN PROCEDURE ###############################
;
;
; ***************************************************************************
;+
;NAME:
;           WVContParamFit
;
; Purpose:  perform the fit of water vapor continuum parameters
;           according to the input data. The H2O continuum is
;           parameterized according to 
;           Rosenkranz, Radio Science, 33(4),  919, 1998 and
;                       Radio Science, 34(4), 1025, 1999
;                       ftp://mesa.mit.edu/phil/lbl_rt
;           The only data considered is pure H2O and for N2-H2O mixtures.
;
; Inputs:   datafile            full path/file name
;           controlfile         name of autogenerated arts control file
;           artspath            full path of the user arts directory,
;                               e.g. '~username/arts/'
;           frange              frequency range for fit (in Hz)
;           H2Otag              full arts H2O tag name
;           H2Omodel            H2O tag associated model
;           H2Ouparam           H2Otag associated userparameters for model='user'
;           H2Olineshape        line shape information vector for H2O lines
;                               e.g. [shape, normalization factor, cutoff]
;           debug               if selected turns additional printed
;                               information on (directed to stddev)
;           catname             spectroscopic line catalog file name
;           catformat           spectroscopic line catalog format
;           catfmin             spectroscopic line catalog lower frequency limit
;           catfmax             spectroscopic line catalog upper frequency limit
;
; Output:   ContParamSet  vector of water vapor continuum parameters
;
; Calling example: WVContParamFit, datafile='ContFitData.aa', $
;                    H2Otag='H2O-PWR98', $
;                    H2Omodel='RosenkranzLines', $
;                    artspath='/home/home01/tkuhn/ARTS/arts/'
;
; History:  2002-01-04    Thomas Kuhn, iup Bremen
;-
; **************************************************************************
;
PRO WVContParamFit, datafile=datafile,             $
                    controlfile=controlfile,       $
                    artspath=artspath,             $
                    artsjobpath=artsjobpath,       $
                    frange=frange,                 $
                    H2Otag=H2Otag,                 $
                    H2Omodel=H2Omodel,             $
                    H2Ouparam=H2Ouparam,           $
                    H2Olineshape=H2Olineshape,     $
                    mirrorline=mirrorline,         $
                    catname=catname,               $
                    catformat=catformat,           $
                    catfmin=catfmin,               $
                    catfmax=catfmax,               $
                    absratio_cs=absratio_cs,       $ ; limit of abs_line/abs_total for C_s
                    FFITMIN=FFITMIN,               $ ; lower frequency limit of abs_tot data
                    FFITMAX=FFITMAX,               $ ; upper frequency limit of abs_tot data
                    absratio_cf=absratio_cf,       $ ; limit of abs_line/abs_total for C_f
                    paramoutfile=paramoutfile,     $ ; ascii output file
                    paramlatexfile=paramlatexfile, $ ; LATEX output file
                    airn2ratio=airn2ratio,         $ ; gamma_air/gamma/N2 ratio
                    SMILESCPU=SMILESCPU              ; select smiles computer
;
; --- CLOSE ALL OPEN UNITS -------------------------------------------------
close, /all

IF NOT KEYWORD_SET(paramlatexfile) THEN paramlatexfile='WVContParamFit.tex'
;openw, latexunit, paramlatexfile, /APPEND, ERROR=err, /get_lun
;IF (err NE 0) THEN BEGIN
;    print,' WVContParamFit> ERROR: parameter output file >>'+paramlatexfile+$
;          '<< can not be opened!'
;    print,' WVContParamFit> error code for output file:',err
;    print,' WVContParamFit> error message:', !ERR_STRING
;    FREE_LUN, latexunit
;ENDIF
;printf, latexunit, $
;  'cat. & cutoff & ls & C^o_s & & x_s & & C^o_f & & x_f & \\ '
;FREE_LUN, latexunit

; --- SELECT THE SMILES COMPUTER -------------------------------------------
if not keyword_set(SMILESCPU) then SMILESCPU='smiles5'


; ---- COMMON BLOCK WITH CONSTANT ------------------------------------------
COMMON UNITCONVERSION, TENLOG10_EULER, dBkm2Npm, Npm2dBkm, Hz2GHz, Pa2hPa
TENLOG10_EULER = 4.3429448190325176 ; = 10*log(e)
; conversion from dB/km to Np/m for absorption units:
dBkm2Npm  = 1.00000e-3 / TENLOG10_EULER;
; conversion from dB/km to 1/m for absorption units:
Npm2dBkm  = (1.00000e0 / dBkm2Npm)
; conversion from Hz to GHz frequency units:
Hz2GHz    = 1.0000e-9
; conversion from Pa to hPa frequency units:
Pa2hPa    = 1.0000e-2
;
; --- set arts job directory -----------------------------------------------
IF NOT KEYWORD_SET(artsjobpath) THEN BEGIN 
    artsjobpath = '/home/home01/tkuhn/ARTS/'
ENDIF
artsjobpath = check_backslash(artsjobpath)
;
; --- CHECK ARTS PATH VARIABLE ---------------------------------------------
IF NOT KEYWORD_SET(artspath) THEN BEGIN
    spawn, 'echo $ARTS_PATH', artspath
ENDIF
artspath = check_backslash(artspath)
;
; --- FILL COMMON BLOCK WITH ARTS SPECIFIC TAG INFO ------------------------
; this will be used in the following by the procedure
; "CreateArtsControlFile" in file "CreateArtsControlFile.pro"
list_arts_tag, path=artspath, ilevel=0, ok=ok
if (ok NE 0) then begin
    print,' WVContParamFit> ERROR in finding the appropriate ARTS tags'
    print,'                 terminate here the job'
    goto, ende
endif
;
; --- get all the measured data of total absorption ------------------------ 
print, ' WVContParamFit> === 1 === read absorption data...'
all_data = read_H2O_data_file( datafile ) ; this is a structure!
;
; --- selection of H2O arts tag ----------------------------------------------
IF NOT KEYWORD_SET(H2Otag)    THEN  H2Otag    = 'H2O'
IF NOT KEYWORD_SET(H2Omodel)  THEN  H2Omodel  = ''
IF NOT KEYWORD_SET(H2Ouparam) THEN  H2Ouparam = [0.0, 0.0, 0.0]
;
; --- selection of H2O line shape --------------------------------------------
IF NOT KEYWORD_SET(H2Olineshape) THEN H2Olineshape = ['no_shape', 'no_norm', '-1']
 ; -> H2O tag = special tag
IF KEYWORD_SET(mirrorline) THEN BEGIN
  IF (STRLOWCASE(H2Olineshape[0]) NE 'no_shape') THEN mirrorline='yes' ELSE mirrorline='NO'
ENDIF ELSE BEGIN
    mirrorline='NO'
ENDELSE
;
; --- selection of H2O line catalog ------------------------------------------
IF NOT KEYWORD_SET(catname)   THEN catname   = ''
IF NOT KEYWORD_SET(catformat) THEN catformat = ''
IF NOT KEYWORD_SET(catfmin)   THEN catfmin   = 0.0e0
IF NOT KEYWORD_SET(catfmax)   THEN catfmax   = 0.0e0
;IF (KEYWORD_SET(catname)       AND $
;   (NOT KEYWORD_SET(catformat) OR $
;    NOT KEYWORD_SET(catfmin)   OR $
;    NOT KEYWORD_SET(catfmax))) THEN BEGIN
;    print,'ERROR in line catalog input information - please check it!'
;    print,'catname: '+catname+' catformat: '+catformat+' catfmin=',catfmin,' catfmax=',catfmax
;    STOP
;ENDIF
;
; --- select fit frequency interval ------------------------------------------
IF NOT KEYWORD_SET(FFITMIN) THEN  FFITMIN=1.0e9  ; [Hz]
IF NOT KEYWORD_SET(FFITMAX) THEN  FFITMAX=1.0e12 ; [Hz]
;
; --- select arts control file name ------------------------------------------
IF NOT KEYWORD_SET(controlfile) THEN BEGIN
    controlfilenamecore ='WVContParamFit'+'_'+H2Otag
ENDIF ELSE BEGIN
    controlfilenamecore = controlfile+'_'+H2Otag
ENDELSE
controlfile = controlfilenamecore+'.arts'
;
; --- calculate the line absorption for each measurement ---------------------
;FOR i = 0, 50 DO BEGIN
FOR i = 0, N_ELEMENTS(all_data[*].Abs)-1 DO BEGIN
;   --- check frequency range ------------------------------------------------
    IF ( ( all_data[i].f LT FFITMIN ) OR $
         ( all_data[i].f GT FFITMAX )   ) THEN GOTO, no_arts_run
;   --- checks ---------------------------------------------------------------
    IF ( all_data[i].WVName   NE 'H2O') THEN GOTO, no_arts_run
    IF ((all_data[i].BuffName NE 'N2') AND $
        (all_data[i].BuffName NE 'XX')) THEN GOTO, no_arts_run
;   --- H2O-N2 measurement ---------------------------------------------------
    IF ((all_data[i].BuffName EQ 'N2') AND $
        (all_data[i].WVName   EQ 'H2O')) THEN BEGIN
;        goto, no_arts_run
        tags               = H2Otag
        tag_models         = H2Omodel
        tag_userparameters = H2Ouparam
        tag_up_n           = [N_ELEMENTS(H2Ouparam)]
        catname            = catname
        catformat          = catformat
        catfmin            = catfmin
        catfmax            = catfmax
        lineshapes         = strarr(3, N_ELEMENTS(tags))
        lineshapes[0:2,0]  = H2Olineshape
        frangemin          = all_data[i].f
        frangemax          = all_data[i].f
        frangesteps        = 2
        ptot               = all_data[i].pwv + all_data[i].pd
        ptzfile            = write_ptz_file(artsjobpath,         $
                                            controlfilenamecore, $
                                            ptot,                $
                                            all_data[i].T)
        vmrtagnames        = tags
        vmrfilenames       = write_vmr_file(artsjobpath,         $
                                            controlfilenamecore, $
                                            tags,                $
                                            ptot,                $
                                            all_data[i].pwv)
        vmrbasename        = controlfile+'_vmr_'
        prangemin          = ptot
        prangemax          = ptot
        prangesteps=2
    ENDIF
;   --- pure H2O measurement -------------------------------------------------
    IF ((all_data[i].BuffName EQ 'XX') AND $
        (all_data[i].WVName   EQ 'H2O')) THEN BEGIN
;        goto, no_arts_run
        tags               = H2Otag
        tag_models         = H2Omodel
        tag_userparameters = H2Ouparam
        tag_up_n           = [N_ELEMENTS(H2Ouparam)]
        catname            = catname
        catformat          = catformat
        catfmin            = catfmin
        catfmax            = catfmax
        lineshapes         = strarr(3, N_ELEMENTS(tags))
        lineshapes[0:2,0]  = H2Olineshape
        frangemin          = all_data[i].f
        frangemax          = all_data[i].f
        frangesteps        = 2
        ptot               = all_data[i].pwv
        ptzfile            = write_ptz_file(artsjobpath,         $
                                            controlfilenamecore, $
                                            ptot,                $
                                            all_data[i].T)
        vmrtagnames        = tags
        vmrfilenames       = write_vmr_file(artsjobpath,         $
                                            controlfilenamecore, $
                                            tags,                $
                                            ptot,                $
                                            all_data[i].pwv)
        vmrbasename        = controlfile+'_vmr_'
        prangemin          = ptot
        prangemax          = ptot
        prangesteps        = 2
    ENDIF
;   --- built arts control file for H2O line absorption calc. ----------------
    print, ' WVContParamFit> === 2 === built arts control file...'
    CreateArtsControlFile, flag=flag, debug=0, $
                          artsjobpath=artsjobpath,  $
                          controlfile=controlfile,  $
;                         -------------------------------------------------
                          tags=tags,                $
                          tag_models=tag_models,    $
                          tag_userparameters=tag_userparameters, $
                          tag_up_n=tag_up_n,                     $
;                         -------------------------------------------------
                          catname=catname,          $
                          catformat=catformat,      $
                          catfmin=catfmin,          $
                          catfmax=catfmax,          $
;                         -------------------------------------------------
                          lineshapes=lineshapes,    $
                          mirrorline=mirrorline,    $
;                         -------------------------------------------------
                          frangemin=frangemin,      $
                          frangemax=frangemax,      $
                          frangesteps=2,            $
;                         -------------------------------------------------
                          ptzfile=ptzfile,          $
;                         -------------------------------------------------
                          vmrtagnames=vmrtagnames,  $
                          vmrfilenames=vmrfilenames,$ 
                          vmrbasename=vmrbasename,  $
;                         -------------------------------------------------
                          prangemin=prangemin,      $
                          prangemax=prangemax,      $
                          prangesteps=2
;
    IF (flag NE 0) THEN BEGIN
        print, 'WVContParamFit> !!! ERROR: arts control file not correctly built'
        print, 'WVContParamFit> !!! arts control file can not be executed, i=',i
        goto, no_arts_run
    ENDIF
;
;   --- run arts job -------------------------------------------------------
    print, ' WVContParamFit> === 3 === run arts control file...'
    CD, artsjobpath, CURRENT=aiidir
;    spawn, 'myarts '+controlfile
    spawn, './myarts_'+SMILESCPU+' '+controlfile
;
;   --- retrieve the calculated line absorption coefficient ----------------
    print, ' WVContParamFit> === 4 === check arts report file...'
    spawn, 'grep Goodby '+artsjobpath+controlfilenamecore+'.rep', findgoodby
    CD, aiidir
    ;;print, 'findgoodby: >>'+findgoodby+'<<'
    IF (STRPOS(findgoodby[0], 'Goodby') LT 0) THEN BEGIN
        print, 'WVContParamFit> !!! ERROR: arts calculation not successful!'
        print, 'WVContParamFit> !!! jump tp the next calculation, present loop index=',i
;        goto, no_arts_run
        STOP
    ENDIF
;
;   --- retrieve the calculated line absorption coefficient ----------------
    print, ' WVContParamFit> === 5 === get line absorption from arts output file...'
    all_data[i].LineMod = H2Otag+'.'+H2Omodel
    file = controlfilenamecore+'.abs.aa'
    all_data[i].labs = read_lineabs_data(artsjobpath, file)
;
;   --- calculate abs_cont = abs_tot - abs_line ----------------------------
    all_data[i].cabs = all_data[i].abs - all_data[i].labs ;
no_arts_run:
ENDFOR



; --- extract some necessary information ------------------------------------
shortcat = STR_SEP(catname, '/')
if (N_ELEMENTS(shortcat) LT 1) then  p='_' else p=shortcat[N_ELEMENTS(shortcat)-1]
catid = STR_SEP(p, '_')
tgid  = STR_SEP(H2Otag, '-')
if (STRLEN(catid[0]) LT 3) then begin
    titletext = STRUPCASE(tgid[1])
endif else begin
    titletext = STRUPCASE(catid[0])
endelse
if (H2Olineshape[2] EQ '-1') then begin
    titletext= titletext+'/nocutoff'
endif else begin
    titletext= titletext+'/cutoff'
endelse



; --- perform H2O-H2O fit of abs_cont --------------------------------------
print, ' WVContParamFit> === 6 === perform the H2O-H2O continuum parameter fit...'
plotfileformat=4
IF NOT KEYWORD_SET(absratio_cs) THEN absratio_cs = 25.00e-2 ; ratio max. allowed abs_line/abs_tot
plotfilename = 'CS_fit_plot_'+H2Otag+'_'+p+'_'+H2Olineshape[0]+'_'+$
               H2Olineshape[2]+'_'+mirrorline+'_'+$
               STRCOMPRESS(string((FFITMIN*1.00e-9),FORMAT='(I3)'))+'-'+$
               STRCOMPRESS(string((FFITMAX*1.00e-9),FORMAT='(I3)'))
; CSfit = [C_s, sigma_C_s, x_s, sigma_x_s, chi2, cov(a,b), correlation(a,b)]
CSfit = do_PWR98_h2oh2o_cont_fit(all_data,                         $
                                 artsjobpath, controlfilenamecore, $
                                 absratio_cs,                      $
                                 FFITMIN, FFITMAX,                 $
                                 plotfilename, plotfileformat,     $
                                 titletext)
IF (N_ELEMENTS(CSfit) LT 7) THEN BEGIN
    print,' WVContParamFit> WARNING! The H2O-H2O continuum parameter fit was NOT successful!'
    goto, ende
ENDIF
spawn,'gzip -f '+plotfilename+'.*'



; --- perform H2O-N2 fit of abs_cont --------------------------------------
print, ' WVContParamFit> === 7 === perform the H2O-N2 continuum parameter fit...'
IF NOT KEYWORD_SET(absratio_cf) THEN absratio_cf = 25.00e-2 ; ratio max. allowed abs_line/abs_tot
plotfilename = 'CF_fit_plot_'+H2Otag+'_'+p+'_'+H2Olineshape[0]+'_'+$
               H2Olineshape[2]+'_'+mirrorline+'_'+$
               STRCOMPRESS(string((FFITMIN*1.00e-9),FORMAT='(I3)'))+'-'+$
               STRCOMPRESS(string((FFITMAX*1.00e-9),FORMAT='(I3)'))
; CFfit = [C_f, sigma_C_f, x_f, sigma_x_f, chi2, cov(a,b), correlation(a,b)]
CFfit = do_PWR98_h2on2_cont_fit(all_data,                              $
                                CSfit[0], CSfit[1], CSfit[2], CSfit[3],$
                                artsjobpath, controlfilenamecore,      $
                                absratio_cf,                           $
                                FFITMIN, FFITMAX,                      $
                                plotfilename, plotfileformat,          $
                                titletext)
IF (N_ELEMENTS(CFfit) LT 7) THEN BEGIN
    print,' WVContParamFit> WARNING! The H2O-N2 continuum parameter fit was NOT successful!'
    goto, ende
ENDIF
spawn,'gzip -f '+plotfilename+'.*'



; --- write parameter to output file --------------------------------------
;; physical units / constants
IF NOT KEYWORD_SET(airn2ratio) THEN airn2ratio=(1.000/1.080) ;; conversion factor for N2 -> air
MPM2artsUnit = 2.3026e-26 ;; unit conversion form MPM units to arts units
;; ouotput file
IF NOT KEYWORD_SET(paramoutfile) THEN paramoutfile='WVContParamFit.CSF'
openw, oounit, paramoutfile, /APPEND, ERROR=err, /get_lun
IF (err NE 0) THEN BEGIN
    print,' WVContParamFit> ERROR: parameter output file >>'+paramoutfile+$
          '<< can not be opened!'
    print,' WVContParamFit> error code for output file:',err
    print,' WVContParamFit> error message:', !ERR_STRING
    FREE_LUN, oounit
ENDIF
printf, oounit,' ' 
printf, oounit,'***********************************************************************'
printf, oounit,' ' 
printf, oounit, GetDatum()
printf, oounit,' H2Otag='+H2Otag+'| H2Omodel='+H2Omodel+'| mirrorline='+mirrorline
printf, oounit,' l-shape='+H2Olineshape[0]+', norm=',H2Olineshape[1]+', cutoff=',H2Olineshape[2]+' Hz'
printf, oounit,' line catalog='+catname
printf, oounit, FORMAT='(A9,E10.3,A11,E10.3)',$
               ' FFITMIN=',FFITMIN,',  FFITMAX=',FFITMAX
printf, oounit, FORMAT='(A13,F5.2,A14,F5.2)',$
               ' absratio_cs=',absratio_cs,',  absratio_cf=',absratio_cf
printf, oounit, FORMAT='(A4,E12.4,A5,E12.4,A21,F6.3,A5,F6.3)',$
               ' Cs=',CSfit[0],' +/- ',CSfit[1],' dB/km/hPa2/GHz2, xs=',CSfit[2],' +/- ',CSfit[3]
if (N_ELEMENTS(CSfit) GT 7) then $
  printf, oounit, FORMAT='(A4,E12.4,A22,F6.2,A7,F6.3,A5,F6.3)',$
  '2Cs=',CSfit[7],' dB/km/hPa2/GHz2 (+/- ',CSfit[9],'%) x_s=',$
  CSfit[8],' +/- ',CSfit[11]
printf, oounit, FORMAT='(A4,E12.4,A17,F6.3)',$
               ' Cs=',(CSfit[0]*MPM2artsUnit),' 1/m/Pa2/Hz2, xs=',CSfit[2]
printf, oounit, FORMAT='(A6,F9.4,A6,F9.4,A7,F9.4)',$
               ' chi2=',CSfit[4],', Cov=',CSfit[5],', Corr=',CSfit[6]
printf, oounit, FORMAT='(A9,E12.4,A5,E12.4,A21,F6.3,A5,F6.3)',$
               ' N2:  Cf=',CFfit[0],' +/- ',CFfit[1],' dB/km/hPa2/GHz2, xf=',CFfit[2],' +/- ',CFfit[3]
if (N_ELEMENTS(CFfit) GT 7) then begin
    printf, oounit,  FORMAT='(A9,E12.4,A22,F6.2,A7,F6.3,A5,F6.3)',$
      '2N2:  Cf=',CFfit[7],' dB/km/hPa2/GHz2 (+/- ',CFfit[9],'%) x_s=',$
      CFfit[8],' +/- ',CFfit[11]
    printf, oounit, FORMAT='(A9,E12.4,A17,F6.3)',$
      '2N2:  Cf=',(CFfit[0]*MPM2artsUnit),$
      ' 1/m/Pa2/Hz2, xf=',CFfit[2]
endif
printf, oounit, FORMAT='(A9,E12.4,A5,E12.4,A21,F6.3,A5,F6.3)',$
               ' air: Cf=',(CFfit[0]*airn2ratio),' +/- ',(CFfit[1]*airn2ratio),$
               ' dB/km/hPa2/GHz2, xf=',CFfit[2],' +/- ',CFfit[3]
printf, oounit, FORMAT='(A9,E12.4,A17,F6.3)',$
               ' air: Cf=',(CFfit[0]*MPM2artsUnit*airn2ratio),$
               ' 1/m/Pa2/Hz2, xf=',CFfit[2]
printf, oounit, FORMAT='(A6,F9.4,A6,F9.4,A7,F9.4)',$
               ' chi2=',CFfit[4],', Cov=',CFfit[5],', Corr=',CFfit[6]
FREE_LUN, oounit



; --- LaTex file -----------------------------------------------------------
openw, latexunit, paramlatexfile, /APPEND, ERROR=err, /get_lun
IF (err NE 0) THEN BEGIN
    print,' WVContParamFit> ERROR: parameter output file >>'+paramlatexfile+$
          '<< can not be opened!'
    print,' WVContParamFit> error code for output file:',err
    print,' WVContParamFit> error message:', !ERR_STRING
    FREE_LUN, latexunit
ENDIF
ppp = STR_SEP(titletext,'/')
slshape = '-'
if ((RSTRPOS(STRUPCASE(H2Olineshape[0]), 'LORENTZ')   GE 0) AND $
    (RSTRPOS(STRUPCASE(H2Olineshape[1]), 'QUADRATIC') GE 0)) then slshape = 'VVW'
if ((RSTRPOS(STRUPCASE(H2Olineshape[0]), 'LORENTZ') GE 0) AND $
    (RSTRPOS(STRUPCASE(H2Olineshape[1]), 'LINEAR')  GE 0)) then slshape = 'L'
if ((RSTRPOS(STRUPCASE(H2Olineshape[0]), 'VOIGT')   GE 0) AND $
    (RSTRPOS(STRUPCASE(H2Olineshape[1]), 'QUADRATIC') GE 0)) then slshape = 'VV'
if ((RSTRPOS(STRUPCASE(H2Olineshape[0]), 'VOIGT') GE 0) AND $
    (RSTRPOS(STRUPCASE(H2Olineshape[1]), 'LINEAR')  GE 0)) then slshape = 'V'
printf, latexunit, $
  FORMAT='(A10,A3,A8,A3,A3,A3,F8.2,A3,F5.1,A5,F5.2,A3,F5.2,A3,F8.2,A3,F5.1,A5,F5.2,A3,F5.2,A3)',$
  ppp[0],' & ',ppp[1], ' & ',slshape,' & ',                $
  (CSfit[7]*1.0000E+08),' & ',CSfit[9],'\%) &',ABS(CSfit[8]),  ' & ',CSfit[11],' & ',$
  (CFfit[7]*1.0000E+09),' & ',CFfit[9],'\%) &',ABS(CFfit[8]),  ' & ',CFfit[11],' \\'
FREE_LUN, latexunit





; --- end of the procedure -------------------------------------------------
ende:
close, /all
print, ' WVContParamFit> end of procedure'
print,'***********************************************************************'
;
END
; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
