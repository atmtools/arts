#ifndef disort_DISORT_h
#define disort_DISORT_h

#ifdef ENABLE_DISORT

#ifdef __cplusplus
extern "C" {
#endif

typedef Index INTEGER;
typedef Index LOGICAL;
typedef Numeric REAL;
typedef float REAL4;
typedef Numeric REAL8;

#define TRUE_ (1)
#define FALSE_ (0)

/* ~~~~~~~~~~~~ */
/* VERSION 1.2 */
/* ~~~~~~~~~~~~ */
/* Subroutine */ 
//int disort_(integer *nlyr, doublereal *dtauc, doublereal *
//            ssalb, doublereal *pmom, doublereal *temper, doublereal *wvnmlo, 
//            doublereal *wvnmhi, logical *usrtau, integer *ntau, doublereal *utau, 
//            integer *nstr, logical *usrang, integer *numu, doublereal *umu, 
//            integer *nphi, doublereal *phi, integer *ibcnd, doublereal *fbeam, 
//            doublereal *umu0, doublereal *phi0, doublereal *fisot, logical *
//            lamber, doublereal *albedo, doublereal *hl, doublereal *btemp, 
//            doublereal *ttemp, doublereal *temis, logical *deltam, logical *plank,
//            logical *onlyfl, doublereal *accur, logical *prnt, char *header, 
//            integer *maxcly, integer *maxulv, integer *maxumu, integer *maxcmu, 
//            integer *maxphi, doublereal *rfldir, doublereal *rfldn, doublereal *
//            flup, doublereal *dfdt, doublereal *uavg, doublereal *uu, doublereal *
//            u0u, doublereal *albmed, doublereal *trnmed, ftnlen header_len);

int disort_(INTEGER *nlyr, REAL *dtauc, REAL *ssalb,
            REAL *pmom, REAL *temper, REAL *wvnmlo, 
            REAL *wvnmhi, LOGICAL *usrtau, INTEGER *ntau, REAL *utau, 
            INTEGER *nstr, LOGICAL *usrang, INTEGER *numu, REAL *umu, 
            INTEGER *nphi, REAL *phi, INTEGER *ibcnd, REAL *fbeam, 
            REAL *umu0, REAL *phi0, REAL *fisot, REAL *intang,
            LOGICAL *lamber, REAL *albedo, REAL *hl, REAL *btemp, 
            REAL *ttemp, REAL *temis, LOGICAL *deltam, LOGICAL *plank,
            LOGICAL *onlyfl, REAL *accur, LOGICAL *prnt, char *header, 
            INTEGER *maxcly, INTEGER *maxulv, INTEGER *maxumu, INTEGER *maxcmu, 
            INTEGER *maxphi, REAL *rfldir, REAL *rfldn,
            REAL *flup, REAL *dfdt, REAL *uavg, REAL *uu,
            REAL *u0u, REAL *albmed, REAL *trnmed);
            
int qgausn_(INTEGER *nn, REAL *gmu, REAL *gwt);
#ifdef __cplusplus
}
#endif

#endif /* ENABLE_DISORT */

#endif /* disort_DISORT_h */

