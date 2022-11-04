/*
 * Copyright (c) 2010-2018 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/** \file sht_com.c
 * \brief Common spherical harmonic tranforms
 */

/// common sht functions and api calls

/* regular transforms */

void SH_to_spat(shtns_cfg shtns, cplx *Qlm, double *Vr) {
	((pf2l)shtns->ftable[SHT_STD][SHT_TYP_SSY])(shtns, Qlm, Vr, shtns->lmax);
}

void spat_to_SH(shtns_cfg shtns, double *Vr, cplx *Qlm) {
	((pf2l)shtns->ftable[SHT_STD][SHT_TYP_SAN])(shtns, Vr, Qlm, shtns->lmax);
}

void SHsphtor_to_spat(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp) {
	((pf4l)shtns->ftable[SHT_STD][SHT_TYP_VSY])(shtns, Slm, Tlm, Vt, Vp, shtns->lmax);
}

void spat_to_SHsphtor(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm) {
	((pf4l)shtns->ftable[SHT_STD][SHT_TYP_VAN])(shtns, Vt, Vp, Slm, Tlm, shtns->lmax);
}

void SHqst_to_spat(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp) {
	((pf6l)shtns->ftable[SHT_STD][SHT_TYP_3SY])(shtns, Qlm, Slm, Tlm, Vr, Vt, Vp, shtns->lmax);
}

void spat_to_SHqst(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm) {
	((pf6l)shtns->ftable[SHT_STD][SHT_TYP_3AN])(shtns, Vr, Vt, Vp, Qlm, Slm, Tlm, shtns->lmax);
}

void SHtor_to_spat(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp) {
	((pf3l)shtns->ftable[SHT_STD][SHT_TYP_GTO])(shtns, Tlm, Vt, Vp, shtns->lmax);
}

void SHsph_to_spat(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp) {
	((pf3l)shtns->ftable[SHT_STD][SHT_TYP_GSP])(shtns, Slm, Vt, Vp, shtns->lmax);
}


/* variable l-truncated transforms */

void SH_to_spat_l(shtns_cfg shtns, cplx *Qlm, double *Vr, int ltr) {
	((pf2l)shtns->ftable[SHT_STD][SHT_TYP_SSY])(shtns, Qlm, Vr, ltr);
}

void spat_to_SH_l(shtns_cfg shtns, double *Vr, cplx *Qlm, int ltr) {
	((pf2l)shtns->ftable[SHT_STD][SHT_TYP_SAN])(shtns, Vr, Qlm, ltr);
}

void SHsphtor_to_spat_l(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, int ltr) {
	((pf4l)shtns->ftable[SHT_STD][SHT_TYP_VSY])(shtns, Slm, Tlm, Vt, Vp, ltr);
}

void spat_to_SHsphtor_l(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, int ltr) {
	((pf4l)shtns->ftable[SHT_STD][SHT_TYP_VAN])(shtns, Vt, Vp, Slm, Tlm, ltr);
}

void SHqst_to_spat_l(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr) {
	((pf6l)shtns->ftable[SHT_STD][SHT_TYP_3SY])(shtns, Qlm, Slm, Tlm, Vr, Vt, Vp, ltr);
}

void spat_to_SHqst_l(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr) {
	((pf6l)shtns->ftable[SHT_STD][SHT_TYP_3AN])(shtns, Vr, Vt, Vp, Qlm, Slm, Tlm, ltr);
}

void SHtor_to_spat_l(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp, int ltr) {
	((pf3l)shtns->ftable[SHT_STD][SHT_TYP_GTO])(shtns, Tlm, Vt, Vp, ltr);
}

void SHsph_to_spat_l(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp, int ltr) {
	((pf3l)shtns->ftable[SHT_STD][SHT_TYP_GSP])(shtns, Slm, Vt, Vp, ltr);
}


/* fixed m, variable l-truncated legendre transforms (no fft) */

void SH_to_spat_ml(shtns_cfg shtns, int im, cplx *Qlm, cplx *Vr, int ltr) {
	((pf2ml)shtns->ftable[SHT_M][SHT_TYP_SSY])(shtns, im, Qlm, Vr, ltr);
}

void spat_to_SH_ml(shtns_cfg shtns, int im, cplx *Vr, cplx *Qlm, int ltr) {
	((pf2ml)shtns->ftable[SHT_M][SHT_TYP_SAN])(shtns, im, Vr, Qlm, ltr);
}

void SHsphtor_to_spat_ml(shtns_cfg shtns, int im, cplx *Slm, cplx *Tlm, cplx *Vt, cplx *Vp, int ltr) {
	((pf4ml)shtns->ftable[SHT_M][SHT_TYP_VSY])(shtns, im, Slm, Tlm, Vt, Vp, ltr);
}

void spat_to_SHsphtor_ml(shtns_cfg shtns, int im, cplx *Vt, cplx *Vp, cplx *Slm, cplx *Tlm, int ltr) {
	((pf4ml)shtns->ftable[SHT_M][SHT_TYP_VAN])(shtns, im, Vt, Vp, Slm, Tlm, ltr);
}

void SHqst_to_spat_ml(shtns_cfg shtns, int im, cplx *Qlm, cplx *Slm, cplx *Tlm, cplx *Vr, cplx *Vt, cplx *Vp, int ltr) {
	((pf6ml)shtns->ftable[SHT_M][SHT_TYP_3SY])(shtns, im, Qlm, Slm, Tlm, Vr, Vt, Vp, ltr);
}

void spat_to_SHqst_ml(shtns_cfg shtns, int im, cplx *Vr, cplx *Vt, cplx *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr) {
	((pf6ml)shtns->ftable[SHT_M][SHT_TYP_3AN])(shtns, im, Vr, Vt, Vp, Qlm, Slm, Tlm, ltr);
}

void SHtor_to_spat_ml(shtns_cfg shtns, int im, cplx *Tlm, cplx *Vt, cplx *Vp, int ltr) {
	((pf3ml)shtns->ftable[SHT_M][SHT_TYP_GTO])(shtns, im, Tlm, Vt, Vp, ltr);
}

void SHsph_to_spat_ml(shtns_cfg shtns, int im, cplx *Slm, cplx *Vt, cplx *Vp, int ltr) {
	((pf3ml)shtns->ftable[SHT_M][SHT_TYP_GSP])(shtns, im, Slm, Vt, Vp, ltr);
}


/* successive scalar + vector for 3D transform (can be faster than simultaneous transform) */
static
void SHqst_to_spat_2l(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr) {
	((pf2l)shtns->ftable[SHT_STD][SHT_TYP_SSY])(shtns, Qlm, Vr, ltr);
	((pf4l)shtns->ftable[SHT_STD][SHT_TYP_VSY])(shtns, Slm, Tlm, Vt, Vp, ltr);
}

static
void spat_to_SHqst_2l(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr) {
	((pf2l)shtns->ftable[SHT_STD][SHT_TYP_SAN])(shtns, Vr, Qlm, ltr);
	((pf4l)shtns->ftable[SHT_STD][SHT_TYP_VAN])(shtns, Vt, Vp, Slm, Tlm, ltr);
}

static
void SHqst_to_spat_2ml(shtns_cfg shtns, int im, cplx *Qlm, cplx *Slm, cplx *Tlm, cplx *Vr, cplx *Vt, cplx *Vp, int ltr) {
	((pf2ml)shtns->ftable[SHT_M][SHT_TYP_SSY])(shtns, im, Qlm, Vr, ltr);
	((pf4ml)shtns->ftable[SHT_M][SHT_TYP_VSY])(shtns, im, Slm, Tlm, Vt, Vp, ltr);
}

static
void spat_to_SHqst_2ml(shtns_cfg shtns, int im, cplx *Vr, cplx *Vt, cplx *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr) {
	((pf2ml)shtns->ftable[SHT_M][SHT_TYP_SAN])(shtns, im, Vr, Qlm, ltr);
	((pf4ml)shtns->ftable[SHT_M][SHT_TYP_VAN])(shtns, im, Vt, Vp, Slm, Tlm, ltr);
}


#if defined(SHT_F77_API)

/*  Fortran 77 api  */

extern shtns_cfg sht_data;

// Fortran API : Call from fortran without the trailing '_'
///@{

	// regular
/// \ingroup fortapi
void shtns_sh_to_spat_(cplx *Qlm, double *Vr) {
	((pf2l)sht_data->ftable[SHT_STD][SHT_TYP_SSY])(sht_data, Qlm, Vr, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_spat_to_sh_(double *Vr, cplx *Qlm) {
	((pf2l)sht_data->ftable[SHT_STD][SHT_TYP_SAN])(sht_data, Vr, Qlm, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_sh_to_spat_cplx_(cplx *alm, cplx *z) {
    SH_to_spat_cplx(sht_data, alm, z);
}

/// \ingroup fortapi
void shtns_spat_cplx_to_sh_(cplx *z, cplx *alm) {
    spat_cplx_to_SH(sht_data, z, alm);
}

/// \ingroup fortapi
void shtns_sphtor_to_spat_(cplx *Slm, cplx *Tlm, double *Vt, double *Vp) {
	((pf4l)sht_data->ftable[SHT_STD][SHT_TYP_VSY])(sht_data, Slm, Tlm, Vt, Vp, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_spat_to_sphtor_(double *Vt, double *Vp, cplx *Slm, cplx *Tlm) {
	((pf4l)sht_data->ftable[SHT_STD][SHT_TYP_VAN])(sht_data, Vt, Vp, Slm, Tlm, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_qst_to_spat_(cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp) {
	((pf6l)sht_data->ftable[SHT_STD][SHT_TYP_3SY])(sht_data, Qlm, Slm, Tlm, Vr, Vt, Vp, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_spat_to_qst_(double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm) {
	((pf6l)sht_data->ftable[SHT_STD][SHT_TYP_3AN])(sht_data, Vr, Vt, Vp, Qlm, Slm, Tlm, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_sph_to_spat_(cplx *Slm, double *Vt, double *Vp) {
	((pf3l)sht_data->ftable[SHT_STD][SHT_TYP_GSP])(sht_data, Slm, Vt, Vp, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_tor_to_spat_(cplx *Tlm, double *Vt, double *Vp) {
	((pf3l)sht_data->ftable[SHT_STD][SHT_TYP_GTO])(sht_data, Tlm, Vt, Vp, sht_data->lmax);
}


	// variable l-truncation
/// \ingroup fortapi
void shtns_sh_to_spat_l_(cplx *Qlm, double *Vr, int* ltr) {
	((pf2l)sht_data->ftable[SHT_STD][SHT_TYP_SSY])(sht_data, Qlm, Vr, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_sh_l_(double *Vr, cplx *Qlm, int* ltr) {
	((pf2l)sht_data->ftable[SHT_STD][SHT_TYP_SAN])(sht_data, Vr, Qlm, *ltr);
}

/// \ingroup fortapi
void shtns_sphtor_to_spat_l_(cplx *Slm, cplx *Tlm, double *Vt, double *Vp, int* ltr) {
	((pf4l)sht_data->ftable[SHT_STD][SHT_TYP_VSY])(sht_data, Slm, Tlm, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_sphtor_l_(double *Vt, double *Vp, cplx *Slm, cplx *Tlm, int* ltr) {
	((pf4l)sht_data->ftable[SHT_STD][SHT_TYP_VAN])(sht_data, Vt, Vp, Slm, Tlm, *ltr);
}

/// \ingroup fortapi
void shtns_qst_to_spat_l_(cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int* ltr) {
	((pf6l)sht_data->ftable[SHT_STD][SHT_TYP_3SY])(sht_data, Qlm, Slm, Tlm, Vr, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_qst_l_(double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int* ltr) {
	((pf6l)sht_data->ftable[SHT_STD][SHT_TYP_3AN])(sht_data, Vr, Vt, Vp, Qlm, Slm, Tlm, *ltr);
}

/// \ingroup fortapi
void shtns_sph_to_spat_l_(cplx *Slm, double *Vt, double *Vp, int* ltr) {
	((pf3l)sht_data->ftable[SHT_STD][SHT_TYP_GSP])(sht_data, Slm, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_tor_to_spat_l_(cplx *Tlm, double *Vt, double *Vp, int* ltr) {
	((pf3l)sht_data->ftable[SHT_STD][SHT_TYP_GTO])(sht_data, Tlm, Vt, Vp, *ltr);
}


	// fixed m, variable l-truncation (no fft, aka Legendre transform)
/// \ingroup fortapi
void shtns_sh_to_spat_ml_(int *im, cplx *Qlm, cplx *Vr, int *ltr) {
	((pf2ml)sht_data->ftable[SHT_M][SHT_TYP_SSY])(sht_data, *im, Qlm, Vr, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_sh_ml_(int *im, cplx *Vr, cplx *Qlm, int *ltr) {
	((pf2ml)sht_data->ftable[SHT_M][SHT_TYP_SAN])(sht_data, *im, Vr, Qlm, *ltr);
}

/// \ingroup fortapi
void shtns_sphtor_to_spat_ml_(int *im, cplx *Slm, cplx *Tlm, cplx *Vt, cplx *Vp, int *ltr) {
	((pf4ml)sht_data->ftable[SHT_M][SHT_TYP_VSY])(sht_data, *im, Slm, Tlm, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_sphtor_ml_(int *im, cplx *Vt, cplx *Vp, cplx *Slm, cplx *Tlm, int *ltr) {
	((pf4ml)sht_data->ftable[SHT_M][SHT_TYP_VAN])(sht_data, *im, Vt, Vp, Slm, Tlm, *ltr);
}

/// \ingroup fortapi
void shtns_qst_to_spat_ml_(int *im, cplx *Qlm, cplx *Slm, cplx *Tlm, cplx *Vr, cplx *Vt, cplx *Vp, int *ltr) {
	((pf6ml)sht_data->ftable[SHT_M][SHT_TYP_3SY])(sht_data, *im, Qlm, Slm, Tlm, Vr, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_qst_ml_(int *im, cplx *Vr, cplx *Vt, cplx *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int *ltr) {
	((pf6ml)sht_data->ftable[SHT_M][SHT_TYP_3AN])(sht_data, *im, Vr, Vt, Vp, Qlm, Slm, Tlm, *ltr);
}

/// \ingroup fortapi
void shtns_tor_to_spat_ml_(int *im, cplx *Tlm, cplx *Vt, cplx *Vp, int *ltr) {
	((pf3ml)sht_data->ftable[SHT_M][SHT_TYP_GTO])(sht_data, *im, Tlm, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_sph_to_spat_ml_(int *im, cplx *Slm, cplx *Vt, cplx *Vp, int *ltr) {
	((pf3ml)sht_data->ftable[SHT_M][SHT_TYP_GSP])(sht_data, *im, Slm, Vt, Vp, *ltr);
}


///@}
#endif
