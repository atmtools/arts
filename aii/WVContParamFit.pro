; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
; ########################### INTERNAL FUNCTIONS ########################### 
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
FUNCTION write_ptz_file, artsjobpath, controlfile, pwv, pd, T
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
printf, unit, string(2.000*(pwv+pd), FORMAT='(E12.6)')+' '+string(T, FORMAT='(E12.6)')+' 0.000000e+00'
printf, unit, string(0.500*(pwv+pd), FORMAT='(E12.6)')+' '+string(T, FORMAT='(E12.6)')+' 1.000000e+00'
;printf, unit, string((pwv+pd), FORMAT='(E12.6)')+' '+string(T, FORMAT='(E12.6)')+' 0.000000e+00'
;printf, unit, string((pwv+pd), FORMAT='(E12.6)')+' '+string(T, FORMAT='(E12.6)')+' 1.000000e+00'
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
FUNCTION write_vmr_file, artsjobpath, controlfile, tags, p
;
; ---- CHECKS ---------------------------------------------
IF (ABS(N_ELEMENTS(tags)-N_ELEMENTS(p)) GT 0) THEN BEGIN
    print, ' write_vmr_file> ERROR in size of tags and pressure detected!'
    print, ' write_vmr_file> pressure and tag vectors must have the same size!'
    print, ' write_vmr_file> tag vector size     :',N_ELEMENTS(tags)
    print, ' write_vmr_file> pressure vector size:',N_ELEMENTS(p)
    FREE_LUN, unit & STOP
ENDIF
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
    ptot = 0.0e0
    FOR k = 0, N_ELEMENTS(p)-1 DO ptot = ptot + p[k]
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
; Viewing contents of file '/net/www/deutsch/idl/idllib/contrib/groupk/regress2.pro'
;
; NAME:
;        REGRESS2
;
; PURPOSE:
;        Multiple linear regression fit.
;        Fit the function:
;        Y(i) = A(0)*X(0,i) + A(1)*X(1,i) + ... +
;             A(Nterms-1)*X(Nterms-1,i)
;
; CATEGORY:
;        G2 - Correlation and regression analysis.
;
; CALLING SEQUENCE:
;        Result = REGRESS(X, Y, W [, YFIT, SIGMA, FTEST, R, RMUL, CHISQ])
;
; INPUTS:
;        X:   array of independent variable data.  X must
;             be dimensioned (Nterms, Npoints) where there are Nterms
;             coefficients to be found (independent variables) and
;             Npoints of samples.
;
;        Y:   vector of dependent variable points, must have Npoints
;             elements.
;
;        W:   vector of weights for each equation, must be a Npoints
;             elements vector.  For instrumental weighting
;             w(i) = 1/standard_deviation(Y(i)), for statistical
;             weighting w(i) = 1./Y(i).   For no weighting set w(i)=1,
;             and also set the RELATIVE_WEIGHT keyword.
;
; OUTPUTS:
;        Function result = coefficients = vector of
;        Nterms elements.  Returned as a column vector.
;
; OPTIONAL OUTPUT PARAMETERS:
;        Yfit:     array of calculated values of Y, Npoints elements.
;
;        Sigma:    Vector of standard deviations for coefficients.
;
;        Ftest:    value of F for test of fit.
;
;        Rmul:     multiple linear correlation coefficient.
;
;        R:        Vector of linear correlation coefficient.
;
;        Chisq:    Reduced chi squared.
;
; KEYWORDS:
;RELATIVE_WEIGHT: if this keyword is non-zero, the input weights
;             (W vector) are assumed to be relative values, and not based
;             on known uncertainties in the Y vector.    This is the case for
;             no weighting W(*) = 1.
;
; PROCEDURE:
;        Adapted from the program REGRES, Page 172, Bevington, Data Reduction
;        and Error Analysis for the Physical Sciences, 1969.
;
; MODIFICATION HISTORY:
;        Written, DMS, RSI, September, 1982.
;        Added RELATIVE_WEIGHT keyword, W. Landsman, August 1991
;        29-AUG-1994:   H.C. Wen - Used simpler, clearer algorithm to determine
;                       fit coefficients. The constant term, A0 is now just one
;                       of the X(iterms,*) vectors, enabling the algorithm to
;                       now return the sigma associated with this constant term.
;                       I also made a special provision for the case when
;                       Nterm = 1; namely, "inverting" a 1x1 matrix, i.e. scalar.
;        26-MAR-1996:   Added the DOUBLE and CHECK keywords to the call to DETERM.
;        02-APR-1996:   Test matrix singularity using second argument in INVERT
;                       instead of call to DETERM.
; see http://www.astro.washington.edu/deutsch-bin/getpro/library27.html?REGRESS2
; **************************************************************************
;
FUNCTION REGRESS2,X,Y,W,Yfit,Sigma,Ftest,R,Rmul,Chisq, RELATIVE_WEIGHT=relative_weight

         On_error,2              ;Return to caller if an error occurs

         NP = N_PARAMS()
         if (NP lt 3) or (NP gt 9) then $
              message,'Must be called with 3-9 parameters: '+$
                      'X, Y, W [, Yfit, Sigma, Ftest, R, RMul, Chisq]'

;  Determine the length of these arrays and the number of sources

         SX        = SIZE( X )
         SY        = SIZE( Y )
         nterm     = SX(1)
         npts      = SY(1)

         if (N_ELEMENTS(W) NE SY(1)) OR $
            (SX(0) NE 2) OR (SY(1) NE SX(2)) THEN $
              message, 'Incompatible arrays.'

         WW   = REPLICATE(1.,nterm) # W
         curv = ( X*WW ) # TRANSPOSE( X )
         beta = X # (Y*W)

         if nterm eq 1 then begin
              sigma  = 1./sqrt(curv)
              X_coeff= beta/curv
         endif else begin
              err     = INVERT( curv, status )

              if (status eq 1) then begin
                   print,'det( Curvature matrix )=0 .. Using REGRESS'
                   X1   = X
                   linechk   = X(0,0) - X(0,fix( npts*randomu(seed) ))
                   if linechk eq 0 then begin
                        print,'Cannot determine sigma for CONSTANT'
                        X1  = X1(1:nterm-1,*)
                   endif

                   coeff = REGRESS( X1,Y,W,Yfit,A0, Sigma,Ftest,R,Rmul,Chisq)

                   if linechk eq 0 then begin
                        coeff     = [A0,reform(coeff)]
                        Sigma     = [ 0,reform(Sigma)]
                        R         = [ 0,R]
                   endif
                   return, coeff
              endif else if (status eq 2) then begin
                print,'WARNING -- small pivot element used in matrix inversion.'
                print,'           significant accuracy probably lost.'
              endif

              diag    = indgen( nterm )
              sigma   = sqrt( err( diag,diag ) )
              X_coeff = err # beta
         endelse

         Yfit     = TRANSPOSE(X_coeff # X)

         dof   = npts - nterm > 1
         Chisq = TOTAL( (Y-Yfit)^2.*W )
         Chisq = Chisq/dof

;   To calculate the "test of fit" parameters, we revert back to the original
;   cryptic routine in REGRESS1. Because the constant term (if any) is now
;   included in the X variable, NPAR = NTERM_regress2 = NTERM_regress1 + 1.

         if nterm eq 1 then goto, SKIP

         SW = TOTAL(W)           ;SUM OF WEIGHTS
         YMEAN = TOTAL(Y*W)/SW   ;Y MEAN
         XMEAN = (X * (REPLICATE(1.,NTERM) # W)) # REPLICATE(1./SW,NPTS)
         WMEAN = SW/NPTS
         WW = W/WMEAN
         ;
         NFREE = NPTS-1          ;DEGS OF FREEDOM
         SIGMAY = SQRT(TOTAL(WW * (Y-YMEAN)^2)/NFREE) ;W*(Y(I)-YMEAN)
         XX = X- XMEAN # REPLICATE(1.,NPTS)      ;X(J,I) - XMEAN(I)
         WX = REPLICATE(1.,NTERM) # WW * XX      ;W(I)*(X(J,I)-XMEAN(I))
         SIGMAX = SQRT( XX*WX # REPLICATE(1./NFREE,NPTS)) ;W(I)*(X(J,I)-XM)*(X(K,I)-XM)
         R = WX #(Y - YMEAN) / (SIGMAX * SIGMAY * NFREE)


         WW1 = WX # TRANSPOSE(XX)

         ARRAY = INVERT(WW1/(NFREE * SIGMAX #SIGMAX))
         A     = (R # ARRAY)*(SIGMAY/SIGMAX)         ;GET COEFFICIENTS

         FREEN = NPTS-NTERM > 1                 ;DEGS OF FREEDOM, AT LEAST 1.

         CHISQ = TOTAL(WW*(Y-YFIT)^2)*WMEAN/FREEN ;WEIGHTED CHI SQUARED
         IF KEYWORD_SET(relative_weight) then varnce = chisq $
                                         else varnce = 1./wmean

         RMUL = TOTAL(A*R*SIGMAX/SIGMAY)         ;MULTIPLE LIN REG COEFF
         IF RMUL LT 1. THEN FTEST = RMUL/(NTERM-1)/ ((1.-RMUL)/FREEN) ELSE FTEST=1.E6
         RMUL = SQRT(RMUL)

SKIP:    return, X_coeff

end
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
!P.POSITION = [0.025, 0.1, 0.825, 0.9]
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
                            linetotmaxratio
;
;
; ---- SAVE SETTINGS -------------------------------------------------------
P_ini = !P
;
; ---- GENERAL CONSTANTS ---------------------------------------------------
COMMON UNITCONVERSION, TENLOG10_EULER, dBkm2Npm, Npm2dBkm, Hz2GHz, Pa2hPa
;
Tref = 300.0E0                  ; [K] for Theta calculation
Rlinetotmax = 7.500e-1          ; ratio of line to total absorption
Rlinetotmin = 0.000e0           ; ratio of line to total absorption
datanlim = 5                    ; minimum of data points for fit
sigma_ratio_abs_tot  = 2.000e-1 ; Gaussian error propagation: sigma_abstot = 20.0% abs_tot 
sigma_ratio_abs_line = 1.000e-1 ; Gaussian error propagation: sigma_absl   = 10.0% abs_l
sigma_ratio_f        = 1.000e-3 ; Gaussian error propagation: sigma_f      =  0.1% f
sigma_ratio_pwv      = 1.000e-3 ; Gaussian error propagation: sigma_PH2O   =  0.1% P_H2O
sigma_ratio_T        = 5.000e-3 ; Gaussian error propagation. sigma_T      =  0.5% T
IF (linetotmaxratio GE 7.500e-1) THEN linetotmaxratio = Rlinetotmax 
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
; --- FIT H2O-H2O DATA -----------------------------------------------------
fitway = 1
IF (fitway EQ 1) THEN BEGIN ;  ----------- LINFIT --------------------------
;   The LINFIT function fits the paired data {xi, yi} to the 
;   linear model, y = A + Bx, by minimizing the Chi-square error statistic. 
;   The result is a two-element vector containing the model parameters [A, B]. 
; o perform fit
    coeff = LINFIT(h2oh2o[0:vecsize-1,1],      $
                   h2oh2o[0:vecsize-1,0],      $
                   SDEV=h2oh2o[0:vecsize-1,3], $
                   PROB=probability , $
                   SIGMA=error,       $
                   CHISQ=CHISQ2,      $
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
    print, 'Chi2=',CHI2_CS[0],', prob. of fit =',probability
    print, '------------------------ C_S LINFIT ------------------------'
ENDIF
;
fitway = 2
IF (fitway EQ 2) THEN BEGIN  ;  ------- simul. fit for f and Theta ---------
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
;    CSTf = REGRESS(XF, YF, Weights, yfit, lnCSOTf, SigmaS, $
;                   FtestS, RS, RmulS, ChisqS)
    CSTf = REGRESS2(XF2, YF, Weights, Yfit, SigmaS, $
                    FtestS, RS, RmulS, ChisqS, RELATIVE_WEIGHT=0)
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
; ==========================={MAKE A PLOT OF C_S}===========================
;
; --- use aii_plot_file for writing into plot output file ------------------
if not keyword_set(plotfilename)   then plotfilename='CS_fit_test_plot'
if not keyword_set(plotfileformat) then plotfileformat=2
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
ymin   = MIN(h2oh2o[0:vecsize-1,0])
ymax   = MAX(h2oh2o[0:vecsize-1,0])
plot, h2oh2o[0:vecsize-1,1],                 $  ; ln(Theta)
      h2oh2o[0:vecsize-1,0],                 $  ; ln(a)
  /NORMAL,                                   $
  xrange=[xlnmin, xlnmax],                     $
  title=TeXtoIDL('fit of H_2O-H_2O data (A. Bauer et al.)', font=0),            $
  xcharsize=1.5,         $
  xtitle=TeXtoIDL('-ln(\Theta)     [1]', font=0),             $
  XTICK_GET = xticks,                        $
  yrange=[ymin, ymax],                       $
  ycharsize=1.5,         $
  ytitle=TeXtoIDL('ln(\alpha / [dB/km/GHz^2/hPa^2])   [1]', font=0),  $
  yTICK_GET = yticks,                        $
  color=colors[0],                           $
  psym=aii_plotsymbols(0),                   $
  symsize=1.25,                              $
  xstyle=2,                                  $
  ystyle=2,                                  $
  /nodata
; T [K] (upper x-axis)
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
FOR i = 0,vecsize-1 DO BEGIN
    yi = ALOG(CS[0,0]) + CS[1,0]*h2oh2o[i,1] ; ln(Cs) + x_s*ln(Theta)
    plots, h2oh2o[i,1],        $  ; ln(Theta)
           yi,                 $  ; ln(Cs) + x_s*ln(Theta)
           color=colors[0],    $
           linestyle=0,        $
           thick=thick,        $
           /CONTINUE
ENDFOR
;
; --- write date and user info ---------------------------------------------
datum = get_userdate_info(1)
xyouts, plotpos[0], plotpos[1]-0.10, datum, CHARSIZE=0.75, CHARTHICK=1.0, /NORMAL
;
; close plot output file
aii_plot_file, action='end', show='yes', print='no', $
               outdir='/home/home01/tkuhn/ARTS'
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
resvec = dblarr(4) ; = [C_s, sigma_C_s, x_s, sigma_x_s]
resvec[0] = CS[0,0]
resvec[1] = 0.5 * ( SIGMA_CS[0,0] + SIGMA_CS[1,0] )
resvec[2] = CS[1,0]
resvec[3] = SIGMA_CS[2,0]
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
                                  linetotmaxratio
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
Tref                 = 3.000E2  ; [K] for Theta calculation
Rlinetotmax          = 7.500e-1 ; ratio of line to total absorption
Rlinetotmin          = 0.000e0  ; ratio of line to total absorption
datanlim             = 5        ; minimum of data points for fit
sigma_ratio_abs_tot  = 2.000e-1 ; Gaussian error propagation: sigma_abstot = 20.0% abs_tot 
sigma_ratio_abs_line = 1.000e-1 ; Gaussian error propagation: sigma_absl   = 10.0% abs_l
sigma_ratio_f        = 1.000e-3 ; Gaussian error propagation: sigma_f      =  0.1% f
sigma_ratio_pwv      = 1.000e-3 ; Gaussian error propagation: sigma_PH2O   =  0.1% P_H2O
sigma_ratio_T        = 5.000e-3 ; Gaussian error propagation. sigma_T      =  0.5% T
IF (linetotmaxratio GE 7.500e-1) THEN linetotmaxratio = Rlinetotmax 
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
                siga[4]  = sigma_ratio_pwv      * (Pa2hPa*data[i].pd)
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
ENDFOR
print,' do_PWR98_cont_fit> # of H2O-N2 data points omitted=',no
IF (j LT datanlim) THEN BEGIN
    print,' do_PWR98_cont_fit> too few data points (#=',j,') for H2O-N2 fit.'
    print,' do_PWR98_cont_fit> jump to the end of this function!'
    goto, read_lineabs_data_ende
ENDIF
print,' do_PWR98_cont_fit> no. of H2O-N2 data points for fit=',j
h2on22 = h2on2[0:j-1, 0:N_ELEMENTS(h2on2[0,*])-1] ; resize array
index   = SORT(h2on22[*,1])                       ; sort in increasing ln(Theta)
vecsize = N_ELEMENTS(h2on22[*,0])                 ; just for simplicity
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
; --- FIT H2O-N2 DATA ------------------------------------------------------
fitway = 1
IF (fitway EQ 1) THEN BEGIN ;  ----------- LINFIT --------------------------
;   The LINFIT function fits the paired data {xi, yi} to the 
;   linear model, y = A + Bx, by minimizing the Chi-square error statistic. 
;   The result is a two-element vector containing the model parameters [A, B]. 
; o perform fit
    coeff = LINFIT(h2on2[0:vecsize-1,1],      $
                   h2on2[0:vecsize-1,0],      $
                   SDEV=h2on2[0:vecsize-1,3], $
                   PROB=probability , $
                   SIGMA=error,       $
                   CHISQ=CHISQ2,      $
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
    print, FORMAT='(A4,E10.3,A1,E10.3,A2,E10.3)','Cf =',CF[0,0],'+',$
           SIGMA_CF[0,0],'/-',SIGMA_CF[0,0]
    print, FORMAT='(A4,F13.6,A3,F13.6)','xf =',CF[1,0],'+/-',SIGMA_CF[2,0]
    print, 'Chi2=',CHI2_CF[0],', prob. of fit =',probability
    print, '------------------------ C_F LINFIT ------------------------'
ENDIF
;
;fitway = 2
IF (fitway EQ 2) THEN BEGIN  ;  ------- simul. fit for f and Theta ---------
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
;    CFTf = REGRESS(XF, YF, Weights, yfit, lnCSOTf, SigmaS, $
;                   FtestS, RS, RmulS, ChisqS)
    CFTf = REGRESS2(XF2, YF, Weights, Yfit, SigmaS, $
                    FtestS, RS, RmulS, ChisqS, RELATIVE_WEIGHT=0)
    CF[0,2]       = EXP(CFTf[0])      ; fit values of C_s
    CF[1,2]       = CFTf[1]           ; fit values of x_s
    CF[2,2]       = CFTf[2]           ; fit values of x_s
    CHI2_CF[2]    = ChisqS            ; chi^2 of the fit
    SIGMA_CF[0,2] = CF[0,2]*(EXP( SigmaS[0])-1.0e0) ; sigma_plus of C_s
    SIGMA_CF[1,2] = CF[0,2]*(EXP(-SigmaS[0])-1.0e0) ; sigma_minus of C_s
    SIGMA_CF[2,2] = SigmaS[1]         ; error in x_s
    SIGMA_CF[3,2] = SigmaS[2]         ; error in x_f
    char_fit_fun[0] = 'REGRESS'
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
; --- set plot environment -------------------------------------------------
if not keyword_set(plotfilename)   then plotfilename='CF_fit_test_plot'
if not keyword_set(plotfileformat) then plotfileformat=2
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
xlnmin = MIN(h2on2[0:vecsize-1,1])              ; min{-ln(Theta)}_i
xlnmax = MAX(h2on2[0:vecsize-1,1])              ; max{-ln(Theta)}_i
ymin   = MIN(h2on2[0:vecsize-1,0])
ymax   = MAX(h2on2[0:vecsize-1,0])
plot, h2on2[0:vecsize-1,1],                  $  ; x=-ln(Theta)
  h2on2[0:vecsize-1,0],                      $  ; y= ln(a)
  /NORMAL,                                   $
  xrange=[xlnmin, xlnmax],                   $
  title=TeXtoIDL('fit of H_2O-N_2 data (A. Bauer et al.)', font=0),  $
  xcharsize=1.5,                             $
  xtitle=TeXtoIDL('-ln(\Theta)     [1]', font=0),                    $
  XTICK_GET = xticks,                        $
  yrange=[ymin, ymax],                       $
  ycharsize=1.5,         $
  ytitle=TeXtoIDL('ln(\alpha / [dB/km/GHz^2/hPa^2])   [1]', font=0), $
  yTICK_GET = yticks,                        $
  color=colors[0],                           $
  psym=aii_plotsymbols(0),                   $
  symsize=1.25,                              $
  xstyle=2,                                  $
  ystyle=2,                                  $
  /nodata
; T [K] (upper x-axis)
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
FOR i = 0,vecsize-1 DO BEGIN
    yi = ALOG(CF[0,0]) + CF[1,0]*h2on2[i,1] ; ln(Cs) + x_s*ln(Theta)
    plots, h2on2[i,1],         $  ; ln(Theta)
           yi,                 $  ; ln(Cs) + x_s*ln(Theta)
           color=colors[0],    $
           linestyle=0,        $
           thick=thick,        $
           /CONTINUE
ENDFOR
;
; --- write date and user info ---------------------------------------------
datum = get_userdate_info(1)
xyouts, plotpos[0], plotpos[1]-0.10, datum, CHARSIZE=0.75, CHARTHICK=1.0, /NORMAL
;
; close plot output file
aii_plot_file, action='end', show='yes', print='no', $
               outdir='/home/home01/tkuhn/ARTS'
;
; ---- END OF FUNCTION -------------------------------------
read_lineabs_data_ende:
; set saved settings back
!P = P_ini
resvec = dblarr(4) ; result vector [C_f, sigma_C_f, x_f, sigma_x_f]]
resvec[0] = CF[0,0]
resvec[1] = 0.5 * ( SIGMA_CF[0,0] + SIGMA_CF[1,0] ) ; take mean as sigma
resvec[2] = CF[1,0]
resvec[3] = SIGMA_CF[2,0]
RETURN, resvec
END
;
;
;
;############################# MAIN PROCEDURE ###############################
;
;
; ***************************************************************************
; Name:     WVContParamFit
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
;
; **************************************************************************
;
PRO WVContParamFit, datafile=datafile, $
                    controlfile=controlfile, $
                    artspath=artspath, $
                    frange=frange, $
                    H2Otag=H2Otag, $
                    H2Omodel=H2Omodel, $
                    H2Ouparam=H2Ouparam, $
                    H2Olineshape=H2Olineshape, $
                    catname=catname, $
                    catformat=catformat, $
                    catfmin=catfmin, $
                    catfmax=catfmax, $
                    absratio_cs=absratio_cs, $
                    absratio_cf=absratio_cf
;
; --- CLOSE ALL OPEN UNITS -------------------------------------------------
close, /all
;
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
artsjobpath = '/home/home01/tkuhn/ARTS/'
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
list_arts_tag, path=artspath, ilevel=0
;
; --- get all the measured data of total absorption ------------------------ 
print, ' WVContParamFit> === 1 === read absorption data...'
all_data = read_H2O_data_file(artsjobpath+datafile) ; this is a structure!
;
; --- selection of H2O arts tag ----------------------------------------------
IF NOT KEYWORD_SET(H2Otag)    THEN  H2Otag    = 'H2O'
IF NOT KEYWORD_SET(H2Omodel)  THEN  H2Omodel  = ''
IF NOT KEYWORD_SET(H2Ouparam) THEN  H2Ouparam = [0.0, 0.0, 0.0]
;
; --- selection of H2O line shape --------------------------------------------
IF NOT KEYWORD_SET(H2Olineshape) THEN H2Olineshape = ['no_shape', 'no_norm', '-1'] ; -> H2O tag = special tag
;
; --- selection of H2O line catalog ------------------------------------------
IF NOT KEYWORD_SET(catname)   THEN catname   = ''
IF NOT KEYWORD_SET(catformat) THEN catformat = ''
IF NOT KEYWORD_SET(catfmin)   THEN catfmin   = 0.0e0
IF NOT KEYWORD_SET(catfmax)   THEN catfmax   = 0.0e0
IF (KEYWORD_SET(catname) AND $
   (NOT KEYWORD_SET(catformat) OR $
    NOT KEYWORD_SET(catfmin)   OR $
    NOT KEYWORD_SET(catfmax))) THEN BEGIN
    print,'ERROR in line catalog input information - please check it!'
    STOP
ENDIF
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
;   --- checks ---------------------------------------------------------------
    IF (all_data[i].WVName   NE 'H2O') THEN GOTO, no_arts_run
    IF ((all_data[i].BuffName NE 'N2') AND $
        (all_data[i].BuffName NE 'XX')) THEN GOTO, no_arts_run
;   --- H2O-N2 measurement ---------------------------------------------------
    IF ((all_data[i].BuffName EQ 'N2') AND $
        (all_data[i].WVName   EQ 'H2O')) THEN BEGIN
;        goto, no_arts_run
        tags               = [H2Otag, 'N2']
        tag_models         = H2Omodel
        tag_userparameters = H2Ouparam
        catname            = catname
        catformat          = catformat
        catfmin            = catfmin
        catfmax            = catfmax
        lineshapes         = strarr(3,N_ELEMENTS(tags))
        lineshapes[0:2,0]    = H2Olineshape
        lineshapes[0:2,1]    = ['no_shape', 'no_norm', '-1']
        frangemin          = all_data[i].f
        frangemax          = all_data[i].f
        frangesteps        = 2
        ptzfile            = write_ptz_file(artsjobpath, controlfilenamecore,$
                                            all_data[i].pwv, $
                                            all_data[i].pd,  $
                                            all_data[i].T)
        vmrtagnames        = tags
        vmrfilenames       = write_vmr_file(artsjobpath, controlfilenamecore,$
                                            tags, $
                                            [all_data[i].pwv, all_data[i].pd])
        vmrbasename        = controlfile+'_vmr_'
        prangemin          = (all_data[i].pwv+all_data[i].pd)
        prangemax          = (all_data[i].pwv+all_data[i].pd)
        prangesteps=2
    ENDIF
;   --- pure H2O measurement -------------------------------------------------
    IF ((all_data[i].BuffName EQ 'XX') AND $
        (all_data[i].WVName   EQ 'H2O')) THEN BEGIN
;        goto, no_arts_run
        tags               = H2Otag
        tag_models         = H2Omodel
        tag_userparameters = H2Ouparam
        catname            = catname
        catformat          = catformat
        catfmin            = catfmin
        catfmax            = catfmax
        lineshapes         = strarr(3,1) & lineshapes[0:2,0] = H2Olineshape
        frangemin          = all_data[i].f
        frangemax          = all_data[i].f
        frangesteps        = 2
        ptzfile            = write_ptz_file(artsjobpath, controlfilenamecore,     $
                                            all_data[i].pwv, $
                                            all_data[i].pd,  $
                                            all_data[i].T)
        vmrtagnames        = tags
        vmrfilenames       = write_vmr_file(artsjobpath, controlfilenamecore,$
                                            tags, $
                                            all_data[i].pwv)
        vmrbasename        = controlfile+'_vmr_'
        prangemin          = (all_data[i].pwv+all_data[i].pd)
        prangemax          = (all_data[i].pwv+all_data[i].pd)
        prangesteps        = 2
    ENDIF
;   --- built arts control file for H2O line absorption calc. ----------------
    print, ' WVContParamFit> === 2 === built arts control file...'
    CreateArtsControlFile, flag=flag, debug=0, $
                          artsjobpath=artsjobpath, $
                          controlfile=controlfile, $
;                         -------------------------------------------------
                          tags=tags, $
                          tag_models=tag_models, $
                          tag_userparameters=tag_userparameters, $
;                         -------------------------------------------------
                          catname=catname, $
                          catformat=catformat, $
                          catfmin=catfmin, $
                          catfmax=catfmax, $
;                         -------------------------------------------------
                          lineshapes=lineshapes, $
;                         -------------------------------------------------
                          frangemin=frangemin, $
                          frangemax=frangemax, $
                          frangesteps=2, $
;                         -------------------------------------------------
                          ptzfile=ptzfile, $
;                         -------------------------------------------------
                          vmrtagnames=vmrtagnames, $
                          vmrfilenames=vmrfilenames,$ 
                          vmrbasename=vmrbasename, $
;                         -------------------------------------------------
                          prangemin=prangemin, $
                          prangemax=prangemax, $
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
    spawn, 'myarts '+controlfile
;
;   --- retrieve the calculated line absorption coefficient ----------------
    print, ' WVContParamFit> === 4 === check arts report file...'
    spawn, 'grep Goodby '+artsjobpath+controlfilenamecore+'.rep', findgoodby
    CD, aiidir
    ;;print, 'findgoodby: >>'+findgoodby+'<<'
    IF (STRPOS(findgoodby[0], 'Goodby') LT 0) THEN BEGIN
        print, 'WVContParamFit> !!! ERROR: arts calculation not successful!'
        print, 'WVContParamFit> !!! jump tp the next calculation, present loop index=',i
        goto, no_arts_run
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
;
; --- perform H2O-H2O fit of abs_cont --------------------------------------
print, ' WVContParamFit> === 6 === perform the H2O-H2O continuum parameter fit...'
IF NOT KEYWORD_SET(absratio_cs) THEN absratio_cs = 3.300e-1 ; ratio max. allowed abs_line/abs_tot
CSfit = do_PWR98_h2oh2o_cont_fit(all_data, artsjobpath, controlfilenamecore, absratio_cs)
IF (N_ELEMENTS(CSfit) NE 4) THEN BEGIN
    print,' WVContParamFit> WARNING! The H2O-H2O continuum parameter fit was NOT successful!'
ENDIF
print,' WVContParamFit> Cs=',CSfit[0],'dB/km/hPa2/GHz2, xs=',CSfit[2]
;
; --- perform H2O-N2 fit of abs_cont --------------------------------------
print, ' WVContParamFit> === 7 === perform the H2O-N2 continuum parameter fit...'
IF NOT KEYWORD_SET(absratio_cf) THEN absratio_cf = 3.300e-1 ; ratio max. allowed abs_line/abs_tot
CFfit = do_PWR98_h2on2_cont_fit(all_data,                              $
                                CSfit[0], CSfit[1], CSfit[2], CSfit[3],$
                                artsjobpath, controlfilenamecore,      $
                                absratio_cf)
IF (N_ELEMENTS(CFfit) NE 4) THEN BEGIN
    print,' WVContParamFit> WARNING! The H2O-N2 continuum parameter fit was NOT successful!'
ENDIF
print,' WVContParamFit> Cs=',CSfit[0],' dB/km/hPa2/GHz2, xs=',CSfit[2]
print,' WVContParamFit> Cf=',CFfit[0],'dB/km/hPa2/GHz2,  xf=',CFfit[2]
;
; --- end of the procedure -------------------------------------------------
ende:
close, /all
print, ' WVContParamFit> end of procedure'
;
END
; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
