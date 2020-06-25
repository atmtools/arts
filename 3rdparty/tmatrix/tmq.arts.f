C   New release including the LAPACK matrix inversion procedure.
C   We thank Cory Davis (University of Edinburgh) for pointing
C   out the possibility of replacing the proprietary NAG matrix
C   inversion routine by the public-domain LAPACK equivalent.

C   CALCULATION OF LIGHT SCATTERING BY POLYDISPERSE, RANDOMLY          
C   ORIENTED PARTICLES OF IDENTICAL AXIALLY SYMMETRIC SHAPE      

C   This version of the code uses EXTENDED PRECISION variables
C   and must be used along with the accompanying files tmq.par.f
C   and lpq.f.

C   Last update 08/06/2005 

C   The code has been developed by Michael Mishchenko at the NASA
C   Goddard Institute for Space Studies, New York. This research
C   was funded by the NASA Radiation Sciences Program.

C   The code can be used without limitations in any not-for-
C   profit scientific research.  We only request that in any
C   publication using the code the source of the code be acknowledged
C   and relevant references (see below) be made.

C   This version of the code is applicable to spheroids,
C   Chebyshev particles, and finite circular cylinders.

C   The computational method is based on the Watermsn's T-matrix
C   approach and is described in detail in the following papers:
C
C   1.  M. I. Mishchenko, Light scattering by randomly oriented
C       axially symmetric particles, J. Opt. Soc. Am. A,
C       vol. 8, 871-882 (1991).
C
C   2.  M. I. Mishchenko, Light scattering by size-shape
C       distributions of randomly oriented axially symmetric
C       particles of a size comparable to a wavelength,
C       Appl. Opt., vol. 32, 4652-4666 (1993).
C
C   3.  M. I. Mishchenko and L. D. Travis, T-matrix computations
C       of light scattering by large spheroidal particles,
C       Opt. Commun., vol. 109, 16-21 (1994).
C
C   4.  M. I. Mishchenko, L. D. Travis, and A. Macke, Scattering
C       of light by polydisperse, randomly oriented, finite
C       circular cylinders, Appl. Opt., vol. 35, 4927-4940 (1996).
C
C   5.  D. J. Wielaard, M. I. Mishchenko, A. Macke, and B. E. Carlson,
C       Improved T-matrix computations for large, nonabsorbing and
C       weakly absorbing nonspherical particles and comparison
C       with geometrical optics approximation, Appl. Opt., vol. 36,
C       4305-4313 (1997).
C                                                                      
C   A general review of the T-matrix approach can be found in
C
C   6.  M. I. Mishchenko, L. D. Travis, and D. W. Mackowski,
C       T-matrix computations of light scattering by nonspherical
C       particles: a review, J. Quant. Spectrosc. Radiat.
C       Transfer, vol. 55, 535-575 (1996).
C
C   The following paper provides a detailed user guide to the
C   T-matrix code:
C
C   7.  M. I. Mishchenko and L. D. Travis, Capabilities and
C       limitations of a current FORTRAN implementation of the
C       T-matrix method for randomly oriented, rotationally
C       symmetric scatterers, J. Quant. Spectrosc. Radiat. Transfer,
C       vol. 60, 309-324 (1998).
C
C   These papers are available in the .pdf format at the web site
C
C   http://www.giss.nasa.gov/~crmim/publications/
C
C   or in hardcopy upon request from Michael Mishchenko
C   Please e-mail your request to crmim@giss.nasa.gov.
C
C   A comprehensive book "Scattering, Absorption, and Emission of
C   Light by Small Particles" (Cambridge University Press, Cambridge,
C   2002) is also available in the .pdf format at the web site
C
C   http://www.giss.nasa.gov/~crmim/books.html

C   Analytical averaging over particle orientations (Ref. 1) makes
C   this method the fastest exact technique currently available.
C   The use of an automatic convergence procedure
C   (Ref. 2) makes the code convenient in massive computations.
C   Ref. 4 describes features specific for finite cylinders as
C   particles with sharp rectangular edges.  Ref. 5 describes further
C   numerical improvements.

C   Since this code uses extended precision variables, it is slower    
C   than the equivalent double-precision code.  The double-precision   
C   code is also available, but it cannot compute as large particles   
C   for the same shape, refractive index, and wavelength as this code. 

C   This is the first part of the full T-matrix code.  The second part,
C   lpq.f, is completely independent of the first part. It contains no
C   T-matrix-specific subroutines and can be compiled separately.
C   The second part of the code replaces the previously implemented
C   standard matrix inversion scheme based on Gaussian elimination
C   by a scheme based on the LU factorization technique.  
C   As described in Ref. 5 above, the use of the LU factorization is
C   especially beneficial for nonabsorbing or weakly absorbing particles.
C   In this code we use the LAPACK implementation of the LU factorization
C   scheme. LAPACK stands for Linear Algebra PACKage. The latter is
C   publicly available at the following internet site:
C
C   http://www.netlib.org/lapack/

C   INPUT PARAMETERS:                                                  
C                                                                      
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.NE.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius                      
C      NDISTR specifies the distribution of equivalent-sphere radii   
C      NDISTR = 1 - modified gamma distribution                        
C           [Eq. (40) of Ref. 7]                                       
C               AXI=alpha                                              
C               B=r_c                                                  
C               GAM=gamma                                              
C      NDISTR = 2 - log-normal distribution                            
C           [Eq. 41) of Ref. 7]                                        
C               AXI=r_g                                                
C               B=[ln(sigma_g)]**2                                    
C      NDISTR = 3 - power law distribution                             
C           [Eq. (42) of Ref. 7]                                       
C                AXI=r_eff (effective radius)                   
C                B=v_eff (effective variance)                
C                Parameters R1 and R2 (see below) are calculated       
C                automatically for given AXI and B
C      NDISTR = 4 - gamma distribution                                 
C           [Eq. (39) of Ref. 7]                                       
C                AXI=a                                                 
C                B=b                                                   
C      NDISTR = 5 - modified power law distribution
C         [Eq. (24) in M. I. Mishchenko et al.,
C         Bidirectional reflectance of flat,
C         optically thick particulate laters: an efficient radiative
C         transfer solution and applications to snow and soil surfaces,
C         J. Quant. Spectrosc. Radiat. Transfer, Vol. 63, 409-432 (1999)].
C                B=alpha
C                                                                      
C      The code computes NPNAX size distributions of the same type     
C      and with the same values of B and GAM in one run.               
C      The parameter AXI varies from AXMAX to AXMAX/NPNAX in steps of  
C      AXMAX/NPNAX.  To compute a single size distribution, use        
C      NPNAX=1 and AXMAX equal to AXI of this size distribution.       
C                                                                      
C      R1 and R2 - minimum and maximum equivalent-sphere
C           radii in the size distribution.
C           They are calculated automatically
C           for the power law distribution with given AXI and B
C           but must be specified for other distributions
C           after the lines
C         
C             DO 600 IAX=1,NPNAX
C                AXI=AXMAX-DAX*DFLOAT(IAX-1)
C         
C           in the main program.
C           For the modified power law distribution (NDISTR=5), the
C           minimum radius is 0, R2 is the maximum radius,
C           and R1 is the intermediate radius at which the
C           n(r)=const dependence is replaced by the power law
C           dependence.
C                                                                      
C      NKMAX.LE.988 is such that NKMAX+2 is the                        
C           number of Gaussian quadrature points used in               
C           integrating over the size distribution for particles with
C           AXI=AXMAX.  For particles with AXI=AXMAX-AXMAX/NPNAX,      
C           AXMAX-2*AXMAX/NPNAX, etc. the number of Gaussian points    
C           linearly decreases.                                       
C           For the modified power law distribution, the number
C           of integration points on the interval [0,R1] is also
C           equal to NKMAX.
C                                                                      
C      LAM - wavelength of light                                       
C      MRR and MRI - real and imaginary parts of the refractive        
C                  index (MRI.GE.0)   
C      EPS and NP - specify the shape of the particles.                
C             For spheroids NP=-1 and EPS is the ratio of the          
C                 horizontal to rotational axes.  EPS is larger than   
C                 1 for oblate spheroids and smaller than 1 for       
C                 prolate spheroids.                                   
C             For cylinders NP=-2 and EPS is the ratio of the          
C                 diameter to the length.                              
C             For Chebyshev particles NP must be positive and 
C                 is the degree of the Chebyshev polynomial, while     
C                 EPS is the deformation parameter                     
C                 [Eq. (33) of Ref. 7].      
C      DDELT - accuracy of the computations                            
C      NPNA - number of equidistant scattering angles (from 0      
C             to 180 deg) for which the scattering matrix is           
C             calculated.                                              
C      NDGS - parameter controlling the number of division points      
C             in computing integrals over the particle surface.        
C             For compact particles, the recommended value is 2.       
C             For highly aspherical particles larger values (3, 4,...) 
C             may be necessary to obtain convergence.                  
C             The code does not check convergence over this parameter. 
C             Therefore, control comparisons of results obtained with  
C             different NDGS-values are recommended.                   
                                                                       

C   OUTPUT PARAMETERS:                                                 
C                                                                      
C      REFF and VEFF - effective radius and effective variance of      
C          the size distribution as defined by Eqs. (43)-(45) of
C          Ref. 7.       
C      CEXT - extinction cross section per particle                    
C      CSCA - scattering cross section per particle                    
C      W - single scattering albedo                                    
C      <cos> - asymmetry parameter of the phase function               
C      ALPHA1,...,BETA2 - coefficients appearing in the expansions
C          of the elements of the scattering matrix in
C          generalized spherical functions
C          [Eqs. (11)-(16) of Ref. 7].
C      F11,...,F44 - elements of the normalized scattering matrix [as      
C          defined by Eqs. (1)-(3) of Ref. 7] versus scattering angle

C   Note that LAM, r_c, r_g, r_eff, a, R1, and R2 must 
C   be given in the same units of length, and that 
C   the dimension of CEXT and CSCA is that of LAM squared (e.g., if        
C   LAM and AXI are given in microns, then CEXT and CSCA are         
C   calculated in square microns).    
                                                                       
C   The physical correctness of the computed results is tested using   
C   the general inequalities derived by van der Mee and Hovenier,      
C   Astron. Astrophys., vol. 228, 559-568 (1990).  Although            
C   the message that the test of van der Mee and Hovenier is satisfied 
C   does not guarantee that the results are absolutely correct,        
C   the message that the test is not satisfied can mean that something 
C   is wrong.                                                          
                                                                       
C   The convergence of the T-matrix method for particles with          
C   different sizes, refractive indices, and aspect ratios can be      
C   dramatically different.  Usually, large sizes and large aspect     
C   ratios cause problems.  The user of this code                      
C   should first experiment with different input parameters in          
C   order to get an idea of the range of applicability of this         
C   technique.  Sometimes decreasing the aspect ratio                  
C   from 3 to 2 can increase the maximum convergent equivalent-        
C   sphere size parameter by a factor of several (Ref. 7).                      
C   The CPU time required rapidly increases with increasing ratio      
C   radius/wavelength and/or with increasing particle asphericity.     
C   This should be taken into account in planning massive computations.
C   Using an optimizing compiler on IBM RISC workstations saves        
C   about 70% of CPU time.                                             
                                                                       
C   Execution can be automatically terminated if dimensions of certain 
C   arrays are not big enough or if the convergence procedure decides  
C   that the accuracy of extended precision variables is insufficient    
C   to obtain a converged T-matrix solution for given particles.       
C   In all cases, a message appears explaining the cause of termination. 
                                                                       
C   The message                                                        
C        "WARNING:  W IS GREATER THAN 1"                               
C   means that the single-scattering albedo exceeds the maximum        
C   possible value 1.  If W is greater than 1 by more than             
C   DDELT, this message can be an indication of numerical              
C   instability caused by extreme values of particle parameters.       
                                                                       
C   The message "WARNING: NGAUSS=NPNG1" means that convergence over    
C   the parameter NG (see Ref. 2) cannot be obtained for the NPNG1     
C   value specified in the PARAMETER statement in the file tmq.par.f. 
C   Often this is not a serious problem, especially for compact         
C   particles.                                                         
                                                                       
C   Larger and/or more aspherical particles may require larger         
C   values of the parameters NPN1, NPN4, and NPNG1 in the file
C   tmq.par.f.  It is recommended to keep NPN1=NPN4+25 and
C   NPNG1=3*NPN1.  Note that the memory requirement increases    
C   as the third power of NPN4. If the memory of                  
C   a computer is too small to accomodate the code in its current    
C   setting, the parameters NPN1, NPN4, and NPNG1 should be
C   decreased. However, this will decrease the maximum size parameter  
C   that can be handled by the code.                                   
                                                                       
C   In some cases any increases of NPN1 will not make the T-matrix     
C   computations convergent.  This means that the particle is just     
C   too "bad" (extreme size parameter and/or extreme aspect ratio      
C   and/or extreme refractive index; see Ref. 7).              
C   The main program contains several PRINT statements which are       
C   currently commentd out.  If uncommented, these statements will     
C   produce numbers which show the convergence rate and can be         
C   used to determine whether T-matrix computations for given particle 
C   parameters will converge at all.   
                                                                       
C   Some of the common blocks are used to save memory rather than      
C   to transfer data.  Therefore, if a compiler produces a warning     
C   message that the lengths of a common block are different in        
C   different subroutines, this is not a real problem.                 
                                                                       
C   The recommended value of DDELT is 0.001.  For bigger values,       
C   false convergence can be obtained.                                 
                                                                       
C   In computations for spheres use EPS=1.000001 instead of EPS=1.     
C   The use of EPS=1 can cause overflows in some rare cases.           
                                                                       
C   To calculate a monodisperse particle, use the options              
C        NPNAX=1                                                       
C        AXMAX=R                                                       
C        B=1D-1   
C        NKMAX=-1                                                      
C        NDISTR=4                                                      
C        ...                                                           
C        DO 600 IAX=1,NPNAX                                            
C           AXI=AXMAX-DAX*DFLOAT(IAX-1)                                
C           R1=0.9999999*AXI                                           
C           R2=1.0000001*AXI                                          
C        ...                                                           
C   where R is the equivalent-sphere radius.                           
                                                                       
C   It is recommended to use the power law rather than the
C   gamma size distribution, because in this case convergent solution
C   can be obtained for larger REFF and VEFF assuming the same
C   maximal R2 (Mishchenko and Travis, Appl. Opt., vol. 33, 7206-7225,
C   1994).
                                                                       
C   For some compilers, DACOS must be raplaced by DARCOS and DASIN     
C   by DARSIN.                                                         
                                                                       
C   If many different size distributions are computed and the
C   refractive index is fixed, then another approach can be more
C   efficient than running this code many times.  Specifically,
C   scattering results should be computed for monodisperse particles
C   with sizes ranging from essentially zero to some maximum value
C   with a small step size (say, 0.02 microns).  These results
C   should be stored on disk and can be used along with spline
C   interpolation to compute scattering by particles with intermediate
C   sizes.  Scattering patterns for monodisperse nonspherical
C   particles in random orientation are (much) smoother than for
C   monodisperse spheres, and spline interpolation usually gives good
C   results. In this way, averaging over any size distribution is a
C   matter of seconds.  For more on size averaging, see Refs. 2 and 4.

C   We would highly appreciate informing me of any problems encountered 
C   with this code.  Please send your message to the following         
C   e-mail address:  CRMIM@GISS.NASA.GOV.                              

C   WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
C   IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
C   INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
C   VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
C   THE AUTHORS AND THEIR ORGANIZATION DISCLAIM ALL LIABILITY FOR
C   ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM. 

      SUBROUTINE TMD(RAT,NDISTR,AXMAX,NPNAX,B,GAM,NKMAX,DEPS,NP,DLAM,
     &     DMRR,
     &     DMRI,DDELT,NPNA,NDGS,R1RAT,R2RAT,QUIET,REFF,VEFF,CEXT,CSCA,
     &     WALB,ASYMM,F11,F22,F33,F44,F12,F34,ERRMSG)

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'tmq.par.f'
      REAL*16 LAM,MRR,MRI,X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2),
     *        AN(NPN1),R(NPNG2),DR(NPNG2),PPI,PIR,PII,P,EPS,A,
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8 XG(1000),WG(1000),TR1(NPN2,NPN2),TI1(NPN2,NPN2),
     &        ALPH1(NPL),ALPH2(NPL),ALPH3(NPL),ALPH4(NPL),BET1(NPL),
     &        BET2(NPL),XG1(2000),WG1(2000),
     &        AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BE1(NPL),BE2(NPL)
      REAL*4
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
      INTEGER*8 QUIET
      CHARACTER ERRMSG*100
      COMMON /CT/ TR1,TI1
      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
      P=QARCOS(-1Q0)
      PIN=P
 
C  OPEN FILES *******************************************************
 
C      OPEN (6,FILE='test')
C      OPEN (10,FILE='tmatr.write')
 
C  INPUT DATA ********************************************************
 
C      RAT=0.5D0
C      NDISTR=3
C      AXMAX=1D0
C      NPNAX=2
C      B=1D-1
C      GAM=0.5D0
C      NKMAX=5
      EPS=DEPS
C      NP=-1
      LAM=DLAM
      MRR=DMRR
      MRI=DMRI
C      DDELT=0.001D0 
C      NPNA=19 
C      NDGS=2
 
      NCHECK=0
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
C     WRITE (6,5454) NCHECK
 5454 FORMAT ('NCHECK=',I1)
      DAX=AXMAX/NPNAX
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (DEPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,DEPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (DEPS,RAT)
C     PRINT 8000, RAT
 8000 FORMAT ('RAT=',F8.6)
      IF(QUIET.EQ.0) THEN
        IF(NP.EQ.-1.AND.EPS.GE.1D0) PRINT 7000,EPS
        IF(NP.EQ.-1.AND.EPS.LT.1D0) PRINT 7001,EPS
        IF(NP.GE.0) PRINT 7100,NP,EPS
        IF(NP.EQ.-2.AND.EPS.GE.1D0) PRINT 7150,EPS
        IF(NP.EQ.-2.AND.EPS.LT.1D0) PRINT 7151,EPS
        PRINT 7400, LAM,MRR,MRI
        PRINT 7200,DDELT
      ENDIF
 7000 FORMAT('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',
     &       I1,'(',F5.2,')')
 7150 FORMAT('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',d8.2)
 7400 FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
      DDELT=0.1D0*DDELT
      DO 600 IAX=1,NPNAX
         AXI=AXMAX-DAX*DFLOAT(IAX-1)
         R1=R1RAT*AXI
         R2=R2RAT*AXI
         NK=MAX(IDINT(AXI*NKMAX/AXMAX+2),1) !MAX call added by Cory to avoid occasional        
         IF (NK.GT.1000) THEN               !problems with nkmax=-1
            WRITE (ERRMSG,8001) NK
            PRINT *,ERRMSG
            RETURN
         ENDIF
C used to be a STOP here
         NK=IDINT(AXI*NKMAX/AXMAX+2)                       
C         IF (NK.GT.1000) PRINT 8001,NK
C         IF (NK.GT.1000) STOP
         IF (NDISTR.EQ.3) CALL POWER (AXI,B,R1,R2)
 8001    FORMAT ('NK=',I4,' I.E., IS GREATER THAN 1000. ',
     &           'EXECUTION TERMINATED.')
         CALL GAUSS (NK,0,0,XG,WG)
         Z1=(R2-R1)*0.5D0
         Z2=(R1+R2)*0.5D0
         Z3=R1*0.5D0
         IF (NDISTR.EQ.5) GO TO 3
         DO I=1,NK
            XG1(I)=Z1*XG(I)+Z2
            WG1(I)=WG(I)*Z1
         ENDDO
         GO TO 4
    3    DO I=1,NK
            XG1(I)=Z3*XG(I)+Z3
            WG1(I)=WG(I)*Z3
         ENDDO
         DO I=NK+1,2*NK
            II=I-NK
            XG1(I)=Z1*XG(II)+Z2
            WG1(I)=WG(II)*Z1
         ENDDO
         NK=NK*2
    4    CALL DISTRB (NK,XG1,WG1,NDISTR,AXI,B,GAM,R1,R2,QUIET,
     &                REFF,VEFF,PIN)
         IF (QUIET.EQ.0) THEN
           PRINT 8002,R1,R2
 8002      FORMAT('R1=',F10.6,'   R2=',F10.6)
           IF (DABS(RAT-1D0).LE.1D-6) PRINT 8003, REFF,VEFF
           IF (DABS(RAT-1D0).GT.1D-6) PRINT 8004, REFF,VEFF
 8003      FORMAT('EQUAL-VOLUME-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
 8004      FORMAT('EQUAL-SURFACE-AREA-SPHERE REFF=',F8.4,
     &          '   VEFF=',F7.4)
           PRINT 7250,NK
 7250      FORMAT('NUMBER OF GAUSSIAN QUADRATURE POINTS ',
     &           'IN SIZE AVERAGING =',I4)
         ENDIF
         DO I=1,NPL
            ALPH1(I)=0D0
            ALPH2(I)=0D0
            ALPH3(I)=0D0
            ALPH4(I)=0D0
            BET1(I)=0D0
            BET2(I)=0D0
         ENDDO   
         CSCAT=0D0
         CEXTIN=0D0
         L1MAX=0
         DO 500 INK=1,NK
            I=NK-INK+1
            A=RAT*XG1(I)
            XEV=2D0*PIN*A/DLAM
            IXXX=XEV+4.05D0*XEV**0.333333D0
            INM1=MAX0(4,IXXX)
            IF (INM1.GE.NPN1) THEN
               WRITE (ERRMSG,7333) NPN1
               PRINT *,ERRMSG
               RETURN
            END IF
C            IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  
     &       '.  EXECUTION TERMINATED')
            QEXT1=0D0
            QSCA1=0D0
            DO 50 NMA=INM1,NPN1
               NMAX=NMA
               MMAX=1
               NGAUSS=NMAX*NDGS
               IF (NGAUSS.GT.NPNG1) THEN
                  WRITE (ERRMSG,7340) NGAUSS
                  PRINT *,ERRMSG
                  RETURN
               END IF
c               IF (NGAUSS.GT.NPNG1) STOP
 7340          FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',
     &                '  EXECUTION TERMINATED')
 7334          FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
 7335 FORMAT('                              NMAX1 =', I3,'  DC2=',D8.2,
     &      '  DC1=',D8.2)
               CALL CONST1(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &                   DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                      DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                          +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
               ENDDO   
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
C              PRINT 7334, NMAX,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               NMIN=DFLOAT(NMAX)/2D0+1D0
               DO 10 N=NMIN,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  DQSCA=DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                      +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  DQEXT=(TR1NN+TR1NN1)*DN1
                  DQSCA=DABS(DQSCA/QSCA)
                  DQEXT=DABS(DQEXT/QEXT)
                  NMAX1=N
                  IF (DQSCA.LE.DDELT.AND.DQEXT.LE.DDELT) GO TO 12
   10          CONTINUE
   12          CONTINUE
c              PRINT 7335, NMAX1,DQSCA,DQEXT
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
               IF (NMA.EQ.NPN1) THEN
                  WRITE (ERRMSG,7333) NPN1
                  PRINT *,ERRMSG
                  RETURN
               END IF
c               IF (NMA.EQ.NPN1) STOP
   50       CONTINUE
   55       NNNGGG=NGAUSS+1
            IF (NGAUSS.EQ.NPNG1) PRINT 7336
            MMAX=NMAX1
            DO 150 NGAUS=NNNGGG,NPNG1
               NGAUSS=NGAUS
               NGGG=2*NGAUSS
 7336          FORMAT('WARNING: NGAUSS=NPNG1')
 7337          FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
               CALL CONST1(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &                   DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                      DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO 104 N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                          +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104          CONTINUE
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
c              PRINT 7337, NGGG,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 155
               IF (NGAUS.EQ.NPNG1) PRINT 7336
  150       CONTINUE
  155       CONTINUE
            QSCA=0D0
            QEXT=0D0
            NNM=NMAX*2
            DO 204 N=1,NNM
               QEXT=QEXT+TR1(N,N)
  204       CONTINUE
            IF (NMAX1.GT.NPN4) THEN
               WRITE (ERRMSG,7550) NMAX1
               PRINT *,ERRMSG
            END IF
 7550       FORMAT ('nmax1 = ',I3, ', i.e. greater than NPN4.',
     &              ' Execution terminated')
C            IF (NMAX1.GT.NPN4) STOP              
            DO 213 N2=1,NMAX1
               NN2=N2+NMAX
               DO 213 N1=1,NMAX1
                  NN1=N1+NMAX
                  ZZ1=TR1(N1,N2)
                  RT11(1,N1,N2)=ZZ1
                  ZZ2=TI1(N1,N2)
                  IT11(1,N1,N2)=ZZ2
                  ZZ3=TR1(N1,NN2)
                  RT12(1,N1,N2)=ZZ3
                  ZZ4=TI1(N1,NN2)
                  IT12(1,N1,N2)=ZZ4
                  ZZ5=TR1(NN1,N2)
                  RT21(1,N1,N2)=ZZ5
                  ZZ6=TI1(NN1,N2)
                  IT21(1,N1,N2)=ZZ6
                  ZZ7=TR1(NN1,NN2)
                  RT22(1,N1,N2)=ZZ7
                  ZZ8=TI1(NN1,NN2)
                  IT22(1,N1,N2)=ZZ8
                  QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                 +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
  213       CONTINUE
C           PRINT 7800,0,DABS(QEXT),QSCA,NMAX
            DO 220 M=1,NMAX1
               CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                     DDR,DRR,DRI,NMAX,NCHECK)
               NM=NMAX-M+1
               NM1=NMAX1-M+1
               M1=M+1
               QSC=0D0
               DO 214 N2=1,NM1
                  NN2=N2+M-1
                  N22=N2+NM
                  DO 214 N1=1,NM1
                     NN1=N1+M-1
                     N11=N1+NM
                     ZZ1=TR1(N1,N2)
                     RT11(M1,NN1,NN2)=ZZ1
                     ZZ2=TI1(N1,N2)
                     IT11(M1,NN1,NN2)=ZZ2
                     ZZ3=TR1(N1,N22)
                     RT12(M1,NN1,NN2)=ZZ3
                     ZZ4=TI1(N1,N22)
                     IT12(M1,NN1,NN2)=ZZ4
                     ZZ5=TR1(N11,N2)
                     RT21(M1,NN1,NN2)=ZZ5
                     ZZ6=TI1(N11,N2)
                     IT21(M1,NN1,NN2)=ZZ6
                     ZZ7=TR1(N11,N22)
                     RT22(M1,NN1,NN2)=ZZ7
                     ZZ8=TI1(N11,N22)
                     IT22(M1,NN1,NN2)=ZZ8
                     QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                       +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
  214          CONTINUE
               NNM=2*NM
               QXT=0D0
               DO 215 N=1,NNM
                  QXT=QXT+TR1(N,N)*2D0
  215          CONTINUE
               QSCA=QSCA+QSC
               QEXT=QEXT+QXT
C              PRINT 7800,M,DABS(QXT),QSC,NMAX
 7800          FORMAT(' m=',I3,'  qxt=',D12.6,'  qsc=',D12.6,
     &                '  nmax=',I3)
  220       CONTINUE
            COEFF1=LAM*LAM*0.5D0/P
            CSCA=QSCA*COEFF1
            CEXT=-QEXT*COEFF1
c           PRINT 7880, NMAX,NMAX1
 7880       FORMAT ('nmax=',I3,'   nmax1=',I3)
            CALL GSP (NMAX1,CSCA,LAM,AL1,AL2,AL3,AL4,BE1,BE2,LMAX)
            L1M=LMAX+1
            L1MAX=MAX(L1MAX,L1M)
            WGII=WG1(I)
            WGI=WGII*CSCA
            DO 250 L1=1,L1M
               ALPH1(L1)=ALPH1(L1)+AL1(L1)*WGI
               ALPH2(L1)=ALPH2(L1)+AL2(L1)*WGI
               ALPH3(L1)=ALPH3(L1)+AL3(L1)*WGI
               ALPH4(L1)=ALPH4(L1)+AL4(L1)*WGI
               BET1(L1)=BET1(L1)+BE1(L1)*WGI
               BET2(L1)=BET2(L1)+BE2(L1)*WGI
  250       CONTINUE
            CSCAT=CSCAT+WGI
            CEXTIN=CEXTIN+CEXT*WGII
C           PRINT 6070, I,NMAX,NMAX1,NGAUSS
 6070       FORMAT(4I6)
  500    CONTINUE
         DO 510 L1=1,L1MAX
            ALPH1(L1)=ALPH1(L1)/CSCAT
            ALPH2(L1)=ALPH2(L1)/CSCAT
            ALPH3(L1)=ALPH3(L1)/CSCAT
            ALPH4(L1)=ALPH4(L1)/CSCAT
            BET1(L1)=BET1(L1)/CSCAT
            BET2(L1)=BET2(L1)/CSCAT
  510    CONTINUE
         WALB=CSCAT/CEXTIN
         CALL HOVENR(L1MAX,ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,QUIET)
         ASYMM=ALPH1(2)/3D0
         IF (QUIET.EQ.0) THEN
            PRINT 9100,CEXTIN,CSCAT,WALB,ASYMM
 9100       FORMAT('CEXT=',D12.6,2X,'CSCA=',D12.6,2X,
     &           2X,'W=',D12.6,2X,'<COS>=',D12.6)
            IF (WALB.GT.1D0) PRINT 9111
 9111       FORMAT ('WARNING: W IS GREATER THAN 1')
c         WRITE (10,580) WALB,L1MAX
C         DO L=1,L1MAX
c            WRITE (10,575) ALPH1(L),ALPH2(L),ALPH3(L),ALPH4(L),
C     &                     BET1(L),BET2(L)           
C         ENDDO   
 575        FORMAT(6D14.7)
 580        FORMAT(D14.8,I8)
         ENDIF
         LMAX=L1MAX-1
         CALL MATR (ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,LMAX,NPNA,QUIET,
     &        F11,F22,F33,F44,F12,F34)
  600 CONTINUE
c      ITIME=MCLOCK()
c      TIME=DFLOAT(ITIME)/6000D0
c      PRINT 1001,TIME
c 1001 FORMAT (' time =',F8.2,' min')
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   INPUT PARAMETERS:                                                 *
C                                                                     *
C   NG = 2*NGAUSS - number of quadrature points on the                *
C                   interval  (-1,1). NGAUSS.LE.NPNG1                 *
C   NMAX,MMAX - maximum dimensions of the arrays.  NMAX.LE.NPN1       *
C               MMAX.LE.NPN1                                          *
C   P - pi                                                            *
C                                                                     *
C   OUTPUT PARAMETERS:                                                *
C                                                                     *
C   X,W - points and weights of the quadrature formula                *
C   AN(N) = n*(n+1)                                                   *
C   ANN(N1,N2) = (1/2)*sqrt((2*n1+1)*(2*n2+1)/(n1*(n1+1)*n2*(n2+1)))  *
C   S(I)=1/sin(arccos(x(i)))                                          *
C   SS(I)=S(I)**2                                                     *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE CONST1 (NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
      IMPLICIT REAL*16 (A-H,O-Z)
      INCLUDE 'tmq.par.f'
      REAL*16 X(NPNG2),W(NPNG2),X1(NPNG1),W1(NPNG1),
     *        X2(NPNG1),W2(NPNG1),
     *        S(NPNG2),SS(NPNG2),
     *        AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
 
      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=QFLOAT(NN)
           D=QSQRT(QFLOAT(2*N+1)/QFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5Q0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE
      NG=2*NGAUSS
      IF (NP.EQ.-2) GO  TO 11
      CALL QGAUSS(NG,0,0,X,W)
      GO TO 19
   11 NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-QCOS(QATAN(EPS))
      CALL QGAUSS(NG1,0,0,X1,W1)
      CALL QGAUSS(NG2,0,0,X2,W2)
      DO 12 I=1,NG1
         W(I)=0.5Q0*(XX+1Q0)*W1(I)
         X(I)=0.5Q0*(XX+1Q0)*X1(I)+0.5Q0*(XX-1Q0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5Q0*XX*W2(I)
         X(I+NG1)=-0.5Q0*XX*X2(I)+0.5Q0*XX
   14 CONTINUE
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE
   19 DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1Q0/(1Q0-Y*Y)
           SS(I)=Y
           SS(NG-I+1)=Y
           Y=QSQRT(Y)
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE
      RETURN
      END
 
C***************************************************************
 
      SUBROUTINE QGAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*16 (A-H,P-Z)
      REAL*16 Z(N),W(N)
      A=1Q0
      B=2Q0
      C=3Q0
      IND=MOD(N,2)
      K=N/2+IND
      F=QFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4Q0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6Q0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0Q0
          NITER=0
          CHECK=1Q-32
   10     PB=1Q0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
c         PRINT 5000, CHECK
          CHECK=CHECK*10Q0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(QABS(PB).GT.CHECK*QABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
 5000 format ('QGAUSS DOES NOT CONVERGE, CHECK=',D10.3)
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   INPUT PARAMETERS:                                                 *
C                                                                     *
C   LAM - wavelength of light                                         *
C   MRR and MRI - real and imaginary parts of the refractive index    *
C   A,EPS,NP - specify shape of the particle                          *
C              (see subroutines RSP1, RSP2, and RSP3)                 *
C   NG = NGAUSS*2 - number of gaussian quadrature points on the       *
C                   interval  (-1,1)                                  *
C   X - gaussian division points                                      *
C   P - pi                                                            *
C                                                                     *
C   OUTPUT INFORMATION:                                               *
C                                                                     *
C   PPI = PI**2 , where PI = (2*P)/LAM (wavenumber)                   *
C   PIR = PPI*MRR                                                     *
C   PII = PPI*MRI                                                     *
C   R and DR - see subroutines RSP1, RSP2, and RSP3                   *
C   DDR=1/(PI*SQRT(R))                                                *
C   DRR+I*DRI=DDR/(MRR+I*MRI)                                         *
C   NMAX - dimension of T(m)-matrices                                 *
C   arrays  J,Y,JR,JI,DJ,DY,DJR,DJI are transferred through           *
C         COMMON /CBESS/ - see subroutine BESS                        *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE VARY (LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,
     *                 R,DR,DDR,DRR,DRI,NMAX)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM,
     *        Z(NPNG2),ZR(NPNG2),ZI(NPNG2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),DDR(NPNG2),
     *        DRR(NPNG2),DRI(NPNG2),
     *        DY(NPNG2,NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      CHARACTER*60 ERRMSG
      NG=NGAUSS*2
      IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,NP,R,DR)
      IF (NP.GE.0) CALL RSP2(X,NG,A,EPS,NP,R,DR)
      IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)
      PI=P*2Q0/LAM
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1Q0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0Q0
      DO 10 I=1,NG
           VV=QSQRT(R(I))
           V=VV*PI
           TA=MAX(TA,V)
           VV=1Q0/V
           DDR(I)=VV
           DRR(I)=PRR*VV
           DRI(I)=PRI*VV
           V1=V*MRR
           V2=V*MRI
           Z(I)=V
           ZR(I)=V1
           ZI(I)=V2
   10 CONTINUE
      IF (NMAX.GT.NPN1) THEN
         WRITE (ERRMSG,9000) NMAX,NPN1
         PRINT *, ERRMSG
      END IF
C     IF (NMAX.GT.NPN1) STOP
 9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
      TB=TA*QSQRT(MRR*MRR+MRI*MRI)
      TB=QMAX1(TB,QFLOAT(NMAX))
      NNMAX1=8.0Q0*QSQRT(QMAX1(TA,QFLOAT(NMAX)))+3Q0
      NNMAX2=(TB+4Q0*(TB**0.33333Q0)+8.0Q0*QSQRT(TB))
      NNMAX2=NNMAX2-NMAX+5
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of the functions R(I)=r(y(I))**2 and                  *
C   DR(I)=((d/dy)r(y))/r(y) and horizontal semi-axis A                *
C   for a spheroid specified by the parameters REV (equal-volume-     *
C   sphere radius) and EPS=A/B (ratio of the semi-axes).              *
C   Y(I)=arccos(X(I))                                                 *
C   1.LE.I.LE.NG                                                      *
C   X - arguments                                                     *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE RSP1 (X,NG,NGAUSS,REV,EPS,NP,R,DR)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NG),R(NG),DR(NG)
      A=REV*EPS**(1Q0/3Q0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1Q0
      DO 50 I=1,NGAUSS
          C=X(I)
          CC=C*C
          SS=1Q0-CC
          S=QSQRT(SS)
          RR=1Q0/(SS+EE*CC)
          R(I)=AA*RR
          R(NG-I+1)=R(I)
          DR(I)=RR*C*S*EE1
          DR(NG-I+1)=-DR(I)
   50 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of the functions R(I)=r(y(i))**2 and                  *
C   DR(I)=((d/dy)r(y))/r(y) and parameter R0 for a Chebyshev          *
C   particle specified by the parameters REV (equal-volume-sphere     *
C   radius), EPS, and N.                                              *
C   Y(I)=arccos(X(I))                                                 *
C   1.LE.I.LE.NG                                                      *
C   X - arguments                                                     *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE RSP2 (X,NG,REV,EPS,N,R,DR)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NG),R(NG),DR(NG)
      DNP=QFLOAT(N)
      DN=DNP*DNP
      DN4=DN*4Q0
      EP=EPS*EPS
      A=1Q0+1.5Q0*EP*(DN4-2Q0)/(DN4-1Q0)
      I=(DNP+0.1Q0)*0.5Q0
      I=2*I
      IF (I.EQ.N) A=A-3Q0*EPS*(1Q0+0.25Q0*EP)/
     *              (DN-1Q0)-0.25Q0*EP*EPS/(9Q0*DN-1Q0)
      R0=REV*A**(-1Q0/3Q0)
      DO 50 I=1,NG
         XI=QARCOS(X(I))*DNP
         RI=R0*(1Q0+EPS*QCOS(XI))
         R(I)=RI*RI
         DR(I)=-R0*EPS*DNP*QSIN(XI)/RI
   50 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of the functions R(I)=R(Y(I))**2 and                  *
C   DR(I)=((d/dy)r(y))/r(y)                                           *
C   for a cylinder specified by the parameters REV (equal-volume-     *
C   sphere radius) and EPS=A/H (ratio of radius to semi-length)       *
C   Y(I)=arccos(X(I))                                                 *
C   1.LE.I.LE.NG                                                      *
C   X - arguments                                                     *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NG),R(NG),DR(NG)
      H=REV*( (2Q0/(3Q0*EPS*EPS))**(1Q0/3Q0) )
      A=H*EPS
      DO 50 I=1,NGAUSS
         CO=-X(I)
         SI=QSQRT(1Q0-CO*CO)
         IF (SI/CO.GT.A/H) GO TO 20
         RAD=H/CO
         RTHET=H*SI/(CO*CO)
         GO TO 30
   20    RAD=A/SI
         RTHET=-A*CO/(SI*SI)
   30    R(I)=RAD*RAD
         R(NG-I+1)=R(I)
         DR(I)=-RTHET/RAD
         DR(NG-I+1)=-DR(I)
   50 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of spherical Bessel functions of the first kind       *
C   J(I,N) = j_n(x) and second kind Y(I,N) = y_n(x)                   *
C   of real-valued argument X(I) and first kind JR(I,N)+I*JI(I,N) =   *
C   = j_n(z) of complex argument Z(I)=XR(I)+I*XI(I), as well as       *
C   the functions                                                     *
C                                                                     *
C   DJ(I,N) = (1/x)(d/dx)(x*j_n(x)) ,                                 *
C   DY(I,N) = (1/x)(d/dx)(x*y_n(x)) ,                                 *
C   DJR(I,N) = Re ((1/z)(d/dz)(z*j_n(x)) ,                            *
C   DJI(I,N) = Im ((1/z)(d/dz)(z*j_n(x)) .                            *
C                                                                     *
C   1.LE.N.LE.NMAX                                                    *
C   NMAX.LE.NPN1                                                      *
C   X,XR,XI - arguments                                               *
C   1.LE.I.LE.NG                                                      *
C   Arrays  J,Y,JR,JI,DJ,DY,DJR,DJI are in                            *
C         COMMON /CBESS/                                              *
C   Parameters NNMAX1 and NMAX2 determine computational accuracy      *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE BESS (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NG),XR(NG),XI(NG),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),
     *        AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1),
     *        ADJ(NPN1),ADY(NPN1),ADJR(NPN1),
     *        ADJI(NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
 
      DO 10 I=1,NG
           XX=X(I)
           CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
           CALL RYB(XX,AY,ADY,NMAX)
           YR=XR(I)
           YI=XI(I)
           CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,NNMAX2)
           DO 10 N=1,NMAX
                J(I,N)=AJ(N)
                Y(I,N)=AY(N)
                JR(I,N)=AJR(N)
                JI(I,N)=AJI(N)
                DJ(I,N)=ADJ(N)
                DY(I,N)=ADY(N)
                DJR(I,N)=ADJR(N)
                DJI(I,N)=ADJI(N)
   10 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of spherical Bessel functions of the first kind j     *
C   of real-valued argument x of orders from 1 to NMAX by using       *
C   backward recursion. Parametr NNMAX determines numerical accuracy. *
C   U - function (1/x)(d/dx)(x*j(x))                                  *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 Y(NMAX),U(NMAX),Z(900)
      L=NMAX+NNMAX
      XX=1Q0/X
      Z(L)=1Q0/(QFLOAT(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1Q0/(QFLOAT(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1Q0/(XX-Z(1))
      Y0=Z0*QCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,NMAX
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-QFLOAT(I)*YI*XX
         Y(I)=YI
   10 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of spherical Bessel functions of the second kind y    *
C   of real-valued argument x of orders from 1 to NMAX by using upward*
C   recursion. V - function (1/x)(d/dx)(x*y(x))                       *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE RYB(X,Y,V,NMAX)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 Y(NMAX),V(NMAX)
      C=QCOS(X)
      S=QSIN(X)
      X1=1Q0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3Q0*X3+X1)*C-3Q0*X2*S
      NMAX1=NMAX-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=QFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,NMAX
  10       V(I)=Y(I-1)-QFLOAT(I)*X1*Y(I)
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of spherical Bessel functions of the first kind       *
C   j=JR+I*JI of complex argument x=XR+I*XI of orders from 1 to NMAX  *
C   by using backward recursion. Parametr NNMAX determines numerical  *
C   accuracy. U=UR+I*UI - function (1/x)(d/dx)(x*j(x))                *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      REAL*16 CYR(NPN1),CYI(NPN1),CZR(2500),CZI(2500),
     *        CUR(NPN1),CUI(NPN1)
      L=NMAX+NNMAX
      XRXI=1Q0/(XR*XR+XI*XI)
      CXXR=XR*XRXI
      CXXI=-XI*XRXI 
      QF=1Q0/QFLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         QF=QFLOAT(2*I1+1)
         AR=QF*CXXR-CZR(I1+1)
         AI=QF*CXXI-CZI(I1+1)
         ARI=1Q0/(AR*AR+AI*AI)
         CZR(I1)=AR*ARI
         CZI(I1)=-AI*ARI
    5 CONTINUE
      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1Q0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=QCOS(XR)*QCOSH(XI)
      CI=-QSIN(XR)*QSINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
      CUR(1)=CU1R
      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I
      DO 10 I=2,NMAX
         QI=QFLOAT(I)
         CYI1R=CYR(I-1)
         CYI1I=CYI(I-1)
         CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
         CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
         AR=CYIR*CXXR-CYII*CXXI
         AI=CYII*CXXR+CYIR*CXXI
         CUIR=CYI1R-QI*AR
         CUII=CYI1I-QI*AI
         CYR(I)=CYIR
         CYI(I)=CYII
         CUR(I)=CUIR
         CUI(I)=CUII
         YR(I)=CYIR
         YI(I)=CYII
         UR(I)=CUIR
         UI(I)=CUII
   10 CONTINUE
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   calculation of the T(0) matrix for an axially symmetric particle  *
C                                                                     *
C   Output information:                                               *
C                                                                     *
C   Arrays  TR1 and TI1 (real and imaginary parts of the              *
C   T(0) matrix) are in COMMON /CT/                                   *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*16 R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      REAL*4 PLUS(NPN6*NPN4*NPN4*8)
      COMMON /TMAT/ PLUS,
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1Q0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2Q0
         ELSE
            CONTINUE
      ENDIF
      SI=1Q0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG (X(I1),NMAX,0,DV1,DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR12=0Q0
                AR21=0Q0
                AI12=0Q0
                AI21=0Q0
                GR12=0Q0
                GR21=0Q0
                GI12=0Q0
                GI21=0Q0
                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0Q0) GO TO 205
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
 
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    URI=DR(I)
                    RRI=RR(I)
 
                    F1=RRI*A22
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
  200           CONTINUE
 
  205           AN12=ANN(N1,N2)*FACTOR 
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
  300 CONTINUE
 
      TPIR=PIR
      TPII=PII
      TPPI=PPI
 
      NM=NMAX
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=0Q0
                TQI(K1,KK2)=0Q0
                TRGQR(K1,KK2)=0Q0
                TRGQI(K1,KK2)=0Q0
 
                TQR(KK1,K2)=0Q0
                TQI(KK1,K2)=0Q0
                TRGQR(KK1,K2)=0Q0
                TRGQI(KK1,K2)=0Q0
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
      CALL TT(NMAX,NCHECK)
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   Calculation of the T(M) matrix, M.GE.1, for an axially symmetric  *
C   particle                                                          *
C                                                                     *
C   Input parameters:                                                 *
C                                                                     *
C   M.GE.1                                                            *
C   NG = NGAUSS*2 - number of gaussian division points on the         *
C        interval  (-1,1)                                             *
C   W - quadrature weights                                            *
C   AN,ANN - see subroutine   CONST                                   *
C   S,SS - see subroutine   CONST                                     *
C   ARRAYS  DV1,DV2,DV3,DV4 are in COMMON /DV/ -                      *
C         see subroutine   DVIG                                       *
C   PPI,PIR,PII - see subroutine   VARY                               *
C   R J DR - see subroutines RSP1 and RSP2                            *
C   DDR,DRR,DRI - see subroutine   VARY                               *
C   NMAX - dimension of the T(M) matrix                               *
C   Arrays  J,Y,JR,JI,DJ,DY,DJR,DJI are in                            *
C        COMMON /CBESS/ - see subroutine   BESS                       *
C                                                                     *
C   Output parameters:                                                *
C                                                                     *
C   Arrays  TR1,TI1 (real and imaginary parts of the T(M) matrix)     *
C   are in COMMON /CT/                                                *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE TMATR (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
      REAL*16 R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      REAL*4 PLUS(NPN6*NPN4*NPN4*8)
      COMMON /TMAT/ PLUS,
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
      MM1=M
      QM=QFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1Q0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2Q0
         ELSE
            CONTINUE
      ENDIF
      SI=1Q0
      NM=NMAX+NMAX
      DO 5 N=1,NM
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG (X(I1),NMAX,M,DV1,DV2)
         DO 25 N=1,NMAX
            SI=SIG(N+M)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           WR=W(I)*R(I)
           DS(I)=S(I)*QM*WR
           DSS(I)=SS(I)*QMM
           RR(I)=WR
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0Q0
                AR12=0Q0
                AR21=0Q0
                AR22=0Q0
                AI11=0Q0
                AI12=0Q0
                AI21=0Q0
                AI22=0Q0
                GR11=0Q0
                GR12=0Q0
                GR21=0Q0
                GR22=0Q0
                GI11=0Q0
                GI12=0Q0
                GI21=0Q0
                GI22=0Q0
                SI=SIG(N1+N2)
 
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
                    AA2=A11*DSS(I)+A22
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1
 
                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI
 
                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII
 
                    URI=DR(I)
                    DSI=DS(I)
                    DSSI=DSS(I)
                    RRI=RR(I)
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0Q0) GO TO 150
 
                    E1=DSI*AA1
                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    IF (NCHECK.EQ.1) GO TO 160
 
  150               F1=RRI*AA2
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    IF (NCHECK.EQ.1) GO TO 200
 
  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
  200           CONTINUE
                AN12=ANN(N1,N2)*FACTOR
                R11(N1,N2)=AR11*AN12
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                R22(N1,N2)=AR22*AN12
                I11(N1,N2)=AI11*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                I22(N1,N2)=AI22*AN12
                RG11(N1,N2)=GR11*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                RG22(N1,N2)=GR22*AN12
                IG11(N1,N2)=GI11*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                IG22(N1,N2)=GI22*AN12
 
  300 CONTINUE
      TPIR=PIR
      TPII=PII
      TPPI=PPI
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
 
                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
      CALL TT(NM,NCHECK)
      RETURN
      END
 
C*****************************************************************
c
c     Calculation of the functiONS
c     DV1(n)=dvig(0,m,n,arccos x)
c     and
c     DV2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)
c     1.LE.N.LE.NMAX
c     0.LE.x.LE.1
 
      SUBROUTINE VIG (X,NMAX,M,DV1,DV2)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 DV1(NPN1), DV2(NPN1)
      A=1Q0
      QS=QSQRT(1Q0-X*X)
      QS1=1Q0/QS
      DO 1 N=1,NMAX
         DV1(N)=0Q0
         DV2(N)=0Q0
    1 CONTINUE
      IF (M.NE.0) GO TO 20
      D1=1Q0
      D2=X  
      DO 5 N=1,NMAX  
         QN=QFLOAT(N)
         QN1=QFLOAT(N+1)
         QN2=QFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1 
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
    5 CONTINUE
      RETURN
   20 QMM=QFLOAT(M*M)
      DO 25 I=1,M
         I2=I*2
         A=A*QSQRT(QFLOAT(I2-1)/QFLOAT(I2))*QS
   25 CONTINUE
      D1=0Q0
      D2=A 
      DO 30 N=M,NMAX
         QN=QFLOAT(N)
         QN2=QFLOAT(2*N+1)
         QN1=QFLOAT(N+1)
         QNM=QSQRT(QN*QN-QMM)
         QNM1=QSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
   30 CONTINUE
      RETURN
      END 
 
C**********************************************************************
C                                                                     *
C   Calculation of the matrix    T = - RG(Q) * (Q**(-1))              *
C                                                                     *
C   Input infortmation is in COMMON /CTT/                             *
C   Output information is in COMMON /CT/                              *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE TT(NMAX,NCHECK)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*16 F(NPN2,NPN2),B(NPN2),WORK(NPN2),COND,
     *       QR(NPN2,NPN2),QI(NPN2,NPN2),
     *       RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *       A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      COMPLEX*32 ZQ(NPN2,NPN2),ZW(NPN2)
      INTEGER IPIV(NPN2),IPVT(NPN2)
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
      NDIM=NPN2
      NNMAX=2*NMAX
 
C	Inversion from LAPACK
 
	DO I=1,NNMAX
	   DO J=1,NNMAX
	      ZQ(I,J)=QCMPLX(QR(I,J),QI(I,J))
	   ENDDO
	ENDDO
	INFO=0
        CALL TMZGETRF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
C       IF (INFO.NE.0) WRITE (6,1100) INFO
        CALL TMZGETRI(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
C       IF (INFO.NE.0) WRITE (6,1100) INFO

 1100   FORMAT ('WARNING:  info=', I2)
	DO I=1,NNMAX
	   DO J=1,NNMAX
	      TR=0D0
	      TI=0D0
	      DO K=1,NNMAX
                 ARR=RGQR(I,K)
                 ARI=RGQI(I,K)
                 AR=ZQ(K,J)
                 AI=QIMAG(ZQ(K,J))
                 TR=TR-ARR*AR+ARI*AI
                 TI=TI-ARR*AI-ARI*AR
              ENDDO
	      TR1(I,J)=TR
	      TI1(I,J)=TI
	   ENDDO
	ENDDO
      RETURN
      END
 
C********************************************************************
C                                                                   *
C   Calculation of the expansion coefficients for the (I,Q,U,V)     *
C   representation.                                                 *
C                                                                   *
C   Input parameters:                                               *
C                                                                   *
C      LAM - wavelength of light                                    *
C      CSCA - scattering cross section                              *
C      TR and TI - elements of the t-matrix. Transferred through    *
C                  COMMON /CTM/                                     *
C      NMAX - dimension of T(M) matrices                            *
C                                                                   *
C   Output infortmation:                                            *
C                                                                   *
C      ALF1,...,ALF4,BET1,BET2 - expansion coefficients             *
C      LMAX - number of coefficients minus 1                        *
C                                                                   *
C********************************************************************
 
      SUBROUTINE GSP(NMAX,CSCA,LAM,ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*8 (A-B,D-H,O-Z),COMPLEX*16 (C)
      REAL*16 LAM
      REAL*8 SSIGN(900)
      REAL*8  CSCA,SSI(NPL),SSJ(NPN1),
     &        ALF1(NPL),ALF2(NPL),ALF3(NPL),
     &        ALF4(NPL),BET1(NPL),BET2(NPL),
     &        TR1(NPL1,NPN4),TR2(NPL1,NPN4),
     &        TI1(NPL1,NPN4),TI2(NPL1,NPN4),
     &        G1(NPL1,NPN6),G2(NPL1,NPN6),
     &        AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4),
     &        FR(NPN4,NPN4),FI(NPN4,NPN4),FF(NPN4,NPN4)
      REAL*4 B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4),
     &       B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4),
     &       D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4),
     &       D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4),
     &       D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4),
     &       PLUS1(NPN6*NPN4*NPN4*8)         
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CIM(NPN1)
 
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
      COMMON /CBESS/ B1R,B1I,B2R,B2I    
      COMMON /SS/ SSIGN
      EQUIVALENCE ( PLUS1(1),TR11(1,1,1) )
      EQUIVALENCE (D1(1,1,1),PLUS1(1)),        
     &            (D2(1,1,1),PLUS1(NPL1*NPN4*NPN4+1)),
     &            (D3(1,1,1),PLUS1(NPL1*NPN4*NPN4*2+1)),
     &            (D4(1,1,1),PLUS1(NPL1*NPN4*NPN4*3+1)), 
     &            (D5R(1,1,1),PLUS1(NPL1*NPN4*NPN4*4+1)) 
 
      CALL FACT
      CALL SIGNUM
      LMAX=2*NMAX
      L1MAX=LMAX+1
      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE
      SSI(1)=1D0
      DO 3 I=1,LMAX
         I1=I+1
         SI=DFLOAT(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE
      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1D0/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            FR(J,I)=CCJ
            FI(J,I)=CCJ*CI
            FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1
 
C *****  CALCULATION OF THE ARRAYS B1 AND B2  *****
 
      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2
 
C     PRINT 3300, B1,B2
 3300 FORMAT (' b1 and b2')
      DO 100 N=1,NMAX
 
C *****  CALCULATION OF THE ARRAYS T1 AND T2  *****
 
 
         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1
               L1=NPN6+M
               TT1=TR11(M1,N,NN)
               TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)
               TT4=TR22(M1,N,NN)
               TT5=TI11(M1,N,NN)
               TT6=TI12(M1,N,NN)
               TT7=TI21(M1,N,NN)
               TT8=TI22(M1,N,NN)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               TR1(L1,NN)=T1+T2
               TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4
               TI2(L1,NN)=T3-T4
               IF(M.EQ.0) GO TO 6
               L1=NPN6-M
               T1=TT1-TT2
               T2=TT3-TT4
               T3=TT5-TT6
               T4=TT7-TT8
               TR1(L1,NN)=T1-T2
               TR2(L1,NN)=T1+T2
               TI1(L1,NN)=T3-T4
               TI2(L1,NN)=T3+T4
    6       CONTINUE
   10    CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****
 
         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1
 
C  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****
 
            CALL CCG(N,N1,NMAX,K1,K2,G1)
            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6
               AAR1=0D0
               AAR2=0D0
               AAI1=0D0
               AAI2=0D0
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)
                  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)
                  RI2=TI2(M1,NN)
                  IF(M.EQ.0) GO TO 12
                  M2=NPN6-M
                  RR1=RR1+TR1(M2,NN)*SIG
                  RI1=RI1+TI1(M2,NN)*SIG
                  RR2=RR2+TR2(M2,NN)*SIG
                  RI2=RI2+TI2(M2,NN)*SIG
   12             AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR
   15       CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****
 
            CALL CCG(N,N1,NMAX,K3,K4,G2)
            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6
            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0
               BBI1=0D0
               BBR2=0D0
               BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2
   30       CONTINUE
   40    CONTINUE
  100 CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****
 
C  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****
 
c     PRINT 3301
 3301 FORMAT(' d1, d2, ...')
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0D0
               DD2=0D0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0
               DD4=0D0
               DD5R=0D0
               DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I
  186       CONTINUE
  190    CONTINUE
  200 CONTINUE
 
C  *****  END OF THE CALCULATION OF THE D-ARRAYS *****
 
C  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****
 
C     PRINT 3303
 3303 FORMAT (' g1, g2, ...')
 
      DK=LAM*LAM/(4D0*CSCA*DACOS(-1D0))
      DO 300 L1=1,L1MAX
         G1L=0D0
         G2L=0D0
         G3L=0D0
         G4L=0D0
         G5LR=0D0
         G5LI=0D0
         L=L1-1
         SL=SSI(L1)*DK
         DO 290 N=1,NMAX
            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)
            IF(NNMAX.LT.NNMIN) GO TO 290
            CALL CCG(N,L,NMAX,K1,K2,G1)
            IF(L.GE.2) CALL CCG(N,L,NMAX,K5,K6,G2)
            NL=N+L
            DO 280  NN=NNMIN,NNMAX
               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)
               DM1=0D0
               DM2=0D0
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE
               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI
               IF(L.LT.2) GO TO 280
               DM3=0D0
               DM4=0D0
               DM5R=0D0
               DM5I=0D0
               MMAX=MIN0(N,NN+2)
               MMIN=MAX0(-N,-NN+2)
               M1MAX=NPN6+MMAX
               M1MIN=NPN6+MMIN
               DO 275 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  SSS1=G2(NPN6-M,NNN)
                  DM3=DM3+SSS1*D3(M1,NN,N)
                  DM4=DM4+SSS1*D4(M1,NN,N)
                  DM5R=DM5R+SSS1*D5R(M1,NN,N)
                  DM5I=DM5I+SSS1*D5I(M1,NN,N)
  275          CONTINUE
               G5LR=G5LR-SSS*DM5R
               G5LI=G5LI-SSS*DM5I
               SSS=G2(NPN4,NNN)*FFN
               G3L=G3L+SSS*DM3
               G4L=G4L+SSS*DM4*SI
  280       CONTINUE
  290    CONTINUE
         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL
         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         LMAX=L
         IF(DABS(G1L).LT.1D-6) GO TO 500
  300 CONTINUE
  500 CONTINUE
      RETURN
      END
 
C****************************************************************
 
C   Calculation of the quantities F(N+1)=0.5*ln(n!)
C   0.LE.N.LE.899
 
      SUBROUTINE FACT
      REAL*8 F(900)
      COMMON /FAC/ F
      F(1)=0D0
      F(2)=0D0
      DO 2 I=3,900
         I1=I-1
         F(I)=F(I1)+0.5D0*DLOG(DFLOAT(I1))
    2 CONTINUE
      RETURN
      END
 
C************************************************************
 
C   Calculation of the array SSIGN(N+1)=sign(n)
C   0.LE.N.LE.899
 
      SUBROUTINE SIGNUM
      REAL*8 SSIGN(900)
      COMMON /SS/ SSIGN
      SSIGN(1)=1D0
      DO 2 N=2,899 
         SSIGN(N)=-SSIGN(N-1)
    2 CONTINUE
      RETURN
      END
 
C******************************************************************
C
C   Calculation of Clebsch-Gordan coefficients
C   (n,m:n1,m1/nn,mm)
C   for given n and n1. m1=mm-m, index mm is found from m as
C   mm=m*k1+k2
C
C   Input parameters :  N,N1,NMAX,K1,K2
C                               N.LE.NMAX
C                               N.GE.1
C                               N1.GE.0
C                               N1.LE.N+NMAX
C   Output parameters : GG(M+NPN6,NN+1) - array of the corresponding
C                                       coefficients
C                               /M/.LE.N
C                               /M1/=/M*(K1-1)+K2/.LE.N1
C                               NN.LE.MIN(N+N1,NMAX)
C                               NN.GE.MAX(/MM/,/N-N1/)
C   If K1=1 and K2=0, then 0.LE.M.LE.N
 
      SUBROUTINE CCG(N,N1,NMAX,K1,K2,GG)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 GG(NPL1,NPN6),CD(0:NPN5),CU(0:NPN5)
      IF(NMAX.LE.NPN4.
     &   AND.0.LE.N1.
     &   AND.N1.LE.NMAX+N.
     &   AND.N.GE.1.
     &   AND.N.LE.NMAX) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCG')
    1 NNF=MIN0(N+N1,NMAX)
      MIN=NPN6-N
      MF=NPN6+N
      IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6
      DO 100 MIND=MIN,MF
         M=MIND-NPN6
         MM=M*K1+K2
         M1=MM-M
         IF(IABS(M1).GT.N1) GO TO 90
         NNL=MAX0(IABS(MM),IABS(N-N1))
         IF(NNL.GT.NNF) GO TO 90
         NNU=N+N1
         NNM=(NNU+NNL)*0.5D0
         IF (NNL.EQ.NNU) NNM=NNL
         CALL CCGIN(N,N1,M,MM,C)
         CU(NNL)=C  
         IF (NNL.EQ.NNF) GO TO 50
         C2=0D0
         C1=C
         DO 7 NN=NNL+1,MIN0(NNM,NNF)
            A=DFLOAT((NN+MM)*(NN-MM)*(N1-N+NN))
            A=A*DFLOAT((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
            A=DFLOAT(4*NN*NN)/A
            A=A*DFLOAT((2*NN+1)*(2*NN-1))
            A=DSQRT(A)
            B=0.5D0*DFLOAT(M-M1)
            D=0D0
            IF(NN.EQ.1) GO TO 5
            B=DFLOAT(2*NN*(NN-1))
            B=DFLOAT((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+
     &               MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN-1)*(NN-1))
            D=D*DFLOAT((2*NN-3)*(2*NN-1))
            D=DFLOAT((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
            D=D*DFLOAT((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
            D=DSQRT(D)
    5       C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50
         CALL DIRECT(N,M,N1,M1,NNU,MM,C)
         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0D0
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DFLOAT((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DFLOAT((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DFLOAT(4*(NN+1)*(NN+1))/A
            A=A*DFLOAT((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DFLOAT(2*(NN+2)*(NN+1))
            B=DFLOAT((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1)
     &               +MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN+2)*(NN+2))
            D=D*DFLOAT((2*NN+5)*(2*NN+3))
            D=DFLOAT((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DFLOAT((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE
   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
c           WRITE (6,*) N,M,N1,M1,NN,MM,GG(MIND,NN+1)
    9    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
 
C*********************************************************************
 
      SUBROUTINE DIRECT (N,M,N1,M1,NN,MM,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(900)
      COMMON /FAC/ F
      C=F(2*N+1)+F(2*N1+1)+F(N+N1+M+M1+1)+F(N+N1-M-M1+1)    
      C=C-F(2*(N+N1)+1)-F(N+M+1)-F(N-M+1)-F(N1+M1+1)-F(N1-M1+1)
      C=DEXP(C)
      RETURN
      END
 
C*********************************************************************
C
C   Calculation of the Clebcsh-Gordan coefficients
C   G=(n,m:n1,mm-m/nn,mm)
C   for given n,n1,m,mm, where NN=MAX(/MM/,/N-N1/)
C                               /M/.LE.N
C                               /MM-M/.LE.N1
C                               /MM/.LE.N+N1
 
      SUBROUTINE CCGIN(N,N1,M,MM,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(900),SSIGN(900)
      COMMON /SS/ SSIGN
      COMMON /FAC/ F
      M1=MM-M
      IF(N.GE.IABS(M).
     &   AND.N1.GE.IABS(M1).
     &   AND.IABS(MM).LE.(N+N1)) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCGIN')
    1 IF (IABS(MM).GT.IABS(N-N1)) GO TO 100
      L1=N
      L2=N1
      L3=M
      IF(N1.LE.N) GO TO 50
      K=N
      N=N1
      N1=K
      K=M
      M=M1
      M1=K
   50 N2=N*2
      M2=M*2
      N12=N1*2
      M12=M1*2
      G=SSIGN(N1+M1+1)
     & *DEXP(F(N+M+1)+F(N-M+1)+F(N12+1)+F(N2-N12+2)-F(N2+2)
     &       -F(N1+M1+1)-F(N1-M1+1)-F(N-N1+MM+1)-F(N-N1-MM+1))
      N=L1
      N1=L2
      M=L3
      RETURN
  100 A=1D0
      L1=M
      L2=MM
      IF(MM.GE.0) GO TO 150
      MM=-MM
      M=-M
      M1=-M1
      A=SSIGN(MM+N+N1+1)
  150 G=A*SSIGN(N+M+1)
     &   *DEXP(F(2*MM+2)+F(N+N1-MM+1)+F(N+M+1)+F(N1+M1+1)
     &        -F(N+N1+MM+2)-F(N-N1+MM+1)-F(-N+N1+MM+1)-F(N-M+1)
     &        -F(N1-M1+1))
      M=L1
      MM=L2
      RETURN
      END
 
C*****************************************************************
 
      SUBROUTINE SAREA (D,RAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (D.GE.1) GO TO 10
      E=DSQRT(1D0-D*D)
      R=0.5D0*(D**(2D0/3D0) + D**(-1D0/3D0)*DASIN(E)/E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
   10 E=DSQRT(1D0-1D0/(D*D))
      R=0.25D0*(2D0*D**(2D0/3D0) + D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E))
     &   /E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
      END
 
c****************************************************************
 
      SUBROUTINE SURFCH (N,E,RAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(60),W(60)
      DN=DFLOAT(N)
      E2=E*E
      EN=E*DN
      NG=60
      CALL GAUSS (NG,0,0,X,W)
      S=0D0
      V=0D0
      DO 10 I=1,NG
         XI=X(I)
         DX=DACOS(XI)
         DXN=DN*DX
         DS=DSIN(DX)
         DSN=DSIN(DXN)
         DCN=DCOS(DXN)
         A=1D0+E*DCN
         A2=A*A
         ENS=EN*DSN
         S=S+W(I)*A*DSQRT(A2+ENS*ENS)
         V=V+W(I)*(DS*A+XI*ENS)*DS*A2
   10 CONTINUE
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0/4D0)**(1D0/3D0)
      RAT=RV/RS
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE SAREAC (EPS,RAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      RAT=(1.5D0/EPS)**(1D0/3D0)
      RAT=RAT/DSQRT( (EPS+2D0)/(2D0*EPS) )
      RETURN
      END
 
C********************************************************************
 
C  Computation of R1 and R2 for a power law size distribution with
C  effective radius A and effective variance B
 
      SUBROUTINE POWER (A,B,R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      COMMON AA,BB
      AA=A
      BB=B
      AX=1D-5
      BX=A-1D-5
      R1=ZEROIN (AX,BX,F,0D0)
      R2=(1D0+B)*2D0*A-R1
      RETURN
      END
 
C***********************************************************************
 
      DOUBLE PRECISION FUNCTION ZEROIN (AX,BX,F,TOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      EPS=1D0
   10 EPS=0.5D0*EPS
      TOL1=1D0+EPS
      IF (TOL1.GT.1D0) GO TO 10
   15 A=AX
      B=BX
      FA=F(A)
      FB=F(B)
   20 C=A
      FC=FA
      D=B-A
      E=D
   30 IF (DABS(FC).GE.DABS(FB)) GO TO 40
   35 A=B
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
   40 TOL1=2D0*EPS*DABS(B)+0.5D0*TOL
      XM=0.5D0*(C-B)
      IF (DABS(XM).LE.TOL1) GO TO 90
   44 IF (FB.EQ.0D0) GO TO 90
   45 IF (DABS(E).LT.TOL1) GO TO 70
   46 IF (DABS(FA).LE.DABS(FB)) GO TO 70
   47 IF (A.NE.C) GO TO 50
   48 S=FB/FA
      P=2D0*XM*S
      Q=1D0-S
      GO TO 60
   50 Q=FA/FC
      R=FB/FC
      S=FB/FA
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
      Q=(Q-1D0)*(R-1D0)*(S-1D0)
   60 IF (P.GT.0D0) Q=-Q
      P=DABS(P)
      IF ((2D0*P).GE.(3D0*XM*Q-DABS(TOL1*Q))) GO TO 70
   64 IF (P.GE.DABS(0.5D0*E*Q)) GO TO 70
   65 E=D
      D=P/Q
      GO TO 80
   70 D=XM
      E=D
   80 A=B
      FA=FB
      IF (DABS(D).GT.TOL1) B=B+D
      IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
      FB=F(B)
      IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20
   85 GO TO 30
   90 ZEROIN=B
      RETURN
      END
 
C***********************************************************************
 
      DOUBLE PRECISION FUNCTION F(R1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON A,B
      R2=(1D0+B)*2D0*A-R1
      F=(R2-R1)/DLOG(R2/R1)-A
      RETURN
      END
 
C**********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
C    N - NUMBER OF POINTS                                             *
C    Z - DIVISION POINTS                                              *
C    W - WEIGHTS                                                      *
C**********************************************************************
 
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END
 
C****************************************************************
 
      SUBROUTINE DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,QUIET,REFF,       
     &                   VEFF,PI)                                    
      IMPLICIT REAL*8 (A-H,O-Z)                                     
      INTEGER*8 QUIET
      REAL*8 YY(NNK),WY(NNK)                                 
      IF (NDISTR.EQ.2) GO TO 100                                
      IF (NDISTR.EQ.3) GO TO 200                                  
      IF (NDISTR.EQ.4) GO TO 300                                 
      IF (NDISTR.EQ.5) GO TO 360
      IF (QUIET.EQ.0) PRINT 1001,AA,BB,GAM
 1001 FORMAT('MODIFIED GAMMA DISTRIBUTION, alpha=',F6.4,'  r_c=',
     &  F6.4,'  gamma=',F6.4)                                    
      A2=AA/GAM                                                  
      DB=1D0/BB
      DO 50 I=1,NNK                                                 
         X=YY(I)                                              
         Y=X**AA                                                 
         X=X*DB
         Y=Y*DEXP(-A2*(X**GAM))                                      
         WY(I)=WY(I)*Y                                              
   50 CONTINUE                                                    
      GO TO 400                                                  
  100 IF (QUIET.EQ.0) PRINT 1002,AA,BB                                                          
 1002 FORMAT('LOG-NORMAL DISTRIBUTION, r_g=',F8.4,
     &         '  [ln(sigma_g)]**2=',F6.4)
      DA=1D0/AA                                                           
      DO 150 I=1,NNK                                                     
         X=YY(I)                                                   
         Y=DLOG(X*DA)                                                
         Y=DEXP(-Y*Y*0.5D0/BB)/X                                
         WY(I)=WY(I)*Y                                             
  150 CONTINUE                                                     
      GO TO 400                                                    
  200 IF (QUIET.EQ.0) PRINT 1003
 1003 FORMAT('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')        
      DO 250 I=1,NNK                                               
         X=YY(I)                                                    
         WY(I)=WY(I)/(X*X*X)                                      
  250 CONTINUE                                                     
      GO TO 400                                                    
  300 IF (QUIET.EQ.0) PRINT 1004,AA,BB
 1004 FORMAT ('GAMMA DISTRIBUTION,  a=',F8.4,'  b=',F6.4)
      B2=(1D0-3D0*BB)/BB                                       
      DAB=1D0/(AA*BB)                                              
      DO 350 I=1,NNK                                               
         X=YY(I)                                                 
         X=(X**B2)*DEXP(-X*DAB)                                  
         WY(I)=WY(I)*X                                 
  350 CONTINUE                                                         
      GO TO 400                                                    
  360 IF (QUIET.EQ.0) PRINT 1005,BB
 1005 FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  alpha=',D10.4)
      DO 370 I=1,NNK
         X=YY(I)
         IF (X.LE.R1) WY(I)=WY(I)
         IF (X.GT.R1) WY(I)=WY(I)*(X/R1)**BB
  370 CONTINUE
  400 CONTINUE                                                       
      SUM=0D0
      DO 450 I=1,NNK
         SUM=SUM+WY(I)
  450 CONTINUE
      SUM=1D0/SUM
      DO 500 I=1,NNK
         WY(I)=WY(I)*SUM
  500 CONTINUE
      G=0D0
      DO 550 I=1,NNK
         X=YY(I)
         G=G+X*X*WY(I)
  550 CONTINUE
      REFF=0D0
      DO 600 I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
  600 CONTINUE
      REFF=REFF/G
      VEFF=0D0
      DO 650 I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
  650 CONTINUE
      VEFF=VEFF/(G*REFF*REFF)
      RETURN                                                                    
      END                                                                       
 
C*************************************************************
 
      SUBROUTINE HOVENR(L1,A1,A2,A3,A4,B1,B2,QUIET)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 QUIET
      REAL*8 A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
      DO 100 L=1,L1
         KONTR=1
         LL=L-1
         DL=DFLOAT(LL)*2D0+1D0
         DDL=DL*0.48D0
         AA1=A1(L)
         AA2=A2(L)
         AA3=A3(L)
         AA4=A4(L)
         BB1=B1(L)
         BB2=B2(L)
         IF(LL.GE.1.AND.DABS(AA1).GE.DL) KONTR=2
         IF(DABS(AA2).GE.DL) KONTR=2
         IF(DABS(AA3).GE.DL) KONTR=2
         IF(DABS(AA4).GE.DL) KONTR=2
         IF(DABS(BB1).GE.DDL) KONTR=2
         IF(DABS(BB2).GE.DDL) KONTR=2
         IF(KONTR.EQ.2) PRINT 3000,LL
         C=-0.1D0
         DO 50 I=1,11
            C=C+0.1D0
            CC=C*C
            C1=CC*BB2*BB2
            C2=C*AA4
            C3=C*AA3
            IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KONTR=2
            IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KONTR=2
            IF(KONTR.EQ.2) PRINT 4000,LL,C
   50    CONTINUE
  100 CONTINUE
      IF((KONTR.EQ.1).AND.(QUIET.EQ.0)) PRINT 2000
 2000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS SATISFIED')
 3000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3)
 4000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3,
     & '   A=',D9.2)
      RETURN
      END
 
C****************************************************************
 
C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
C    COEFFICIENTS
 
C    A1,...,B2 - EXPANSION COEFFICIENTS
C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
C    N - NUMBER OF SCATTERING ANGLES
C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MATR(A1,A2,A3,A4,B1,B2,LMAX,NPNA,QUIET,F11,F22,F33,F44,
     &                F12,F34)
      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL),F11(NPNA),
     &     F22(NPNA),F33(NPNA),F44(NPNA),F12(NPNA),F34(NPNA)
      INTEGER QUIET
      N=NPNA
      DN=1D0/DFLOAT(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      L1MAX=LMAX+1
      TB=-DB
      TAA=-DA
      IF (QUIET.EQ.0) THEN
         PRINT 1000
 1000    FORMAT(' ')
         PRINT 1001
 1001    FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',
     &        6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
         DO 10 L1=1,L1MAX
            L=L1-1
            PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
 10      CONTINUE
 1002    FORMAT(' ',I3,6F12.5)
         PRINT 1000
         PRINT 1003
 1003    FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',
     &        8X,'F44',8X,'F12',8X,'F34')
      ENDIF
      D6=DSQRT(6D0)*0.25D0
      DO 500 I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)
         F11(I1)=0D0
         F2=0D0
         F3=0D0
         F44(I1)=0D0
         F12(I1)=0D0
         F34(I1)=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11(I1)=F11(I1)+A1(L1)*PP1
            F44(I1)=F44(I1)+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12(I1)=F12(I1)+B1(L1)*PP4
            F34(I1)=F34(I1)+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE
         F22(I1)=(F2+F3)*0.5D0
         F33(I1)=(F2-F3)*0.5D0
C        F22=F22/F11
C        F33=F33/F11
C        F44=F44/F11
C        F12=-F12/F11
C        F34=F34/F11
         IF (QUIET.EQ.0) PRINT 1004,TB,F11(I1),F22(I1),F33(I1),F44(I1),
     &        F12(I1),F34(I1)
  500 CONTINUE
      IF (QUIET.EQ.0) THEN
         PRINT 1000 
 1004    FORMAT(' ',F6.2,6F11.4)
      ENDIF
      RETURN
      END
