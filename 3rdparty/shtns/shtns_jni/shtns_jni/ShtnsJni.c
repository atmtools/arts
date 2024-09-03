/*
 * Copyright (c) 2017 Julien Pierret
 * 
 * julien.pierret@gmail.com
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

/** 
 A JNI Wrapper to use SHTns from JAVA.
**/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>

#include <shtns.h>
#include <stdint.h>
#include "ShtnsJni.h"

JNIEXPORT jlong JNICALL Java_com_mrsharky_shtns_Shtns_initializeGaussian(JNIEnv *env, jobject thisObj
, jint jlmax, jint jmmax, jint periodicity, jint jnlat, jint jnlong, jint numThreads) {  
    shtns_cfg shtns;
       
    const int lmax = jlmax;             // maximum degree of spherical harmonics  
    const int mmax = jmmax;             // maximum order of spherical harmonics
    const int mres = periodicity;	// periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)
    const int nlat = jnlat;		// number of points in the latitude direction  (constraint: nlat >= lmax+1)
    const int nphi = jnlong;		// number of points in the longitude direction (constraint: nphi >= 2*mmax+1)
    
    shtns_verbose(0);			// displays informations during initialization.
    shtns_use_threads(numThreads);	// enable multi-threaded transforms (if supported).
    shtns = shtns_init( sht_gauss, lmax, mmax, mres, nlat, nphi );
    return (long) shtns;
}

JNIEXPORT jlong JNICALL Java_com_mrsharky_shtns_Shtns_getNlm(JNIEnv *env, jobject thisObj, jlong cfg) {
    shtns_cfg shtns = (shtns_cfg) cfg;
    return shtns->nlm;
}

JNIEXPORT jdoubleArray JNICALL Java_com_mrsharky_shtns_Shtns_spatialToSpectral(JNIEnv *env, jobject thisObj, jlong cfg, jdoubleArray data) {
    shtns_cfg shtns = (shtns_cfg) cfg;
    jdouble *Sh = (*env)->GetDoubleArrayElements(env, data, NULL);
    
    //printf("lmax= %i, mmax=%i, mres=%i, nlat=%i, nphi=%i\n", shtns->lmax, shtns->mmax, shtns->mres, shtns->nlat, shtns->nphi);
    
    complex double *Slm;	// spherical harmonics coefficients (l,m space): complex numbers.
    double *Th;
    //	shtns_set_grid(shtns, sht_gauss, 0.0, nlat, nphi);
    long int NLM = shtns->nlm;
    
    // Memory allocation : the use of shtns_malloc is highly recommended because we need proper alignement.
    // allocate spatial fields.
    Th = (double *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(double));

    // allocate SH representations.
    Slm = (complex double *) shtns_malloc( NLM * sizeof(complex double));

    // go to spectral space: Sh -> Slm
    spat_to_SH(shtns, Sh, Slm);			// note that this destroys the spatial data in Sh.

    int arraySize = shtns->nlm;
    jdouble spectral[2*arraySize];              // Double the size of the return array. We need to store both real & imaginary values in it
    
    for (int j = 0; j < arraySize; j++) {
        spectral[2*j] = creal(Slm[j]);          // Real
        spectral[2*j+1] = cimag(Slm[j]);        // Imaginary
    }
    
    jdoubleArray spectralArray = (*env)->NewDoubleArray(env, 2*arraySize);         // allocate
    (*env)->SetDoubleArrayRegion(env, spectralArray, 0 , 2*arraySize, spectral);   // copy
    
    // Cleanup
    (*env)->ReleaseDoubleArrayElements(env, data, Sh, 0);
    return spectralArray;
}

JNIEXPORT jdoubleArray JNICALL Java_com_mrsharky_shtns_Shtns_spectralToSpatial(JNIEnv *env, jobject thisObj, jlong cfg, jdoubleArray real, jdoubleArray imag) {
    shtns_cfg shtns = (shtns_cfg) cfg;   
    jdouble *Slm_real = (*env)->GetDoubleArrayElements(env, real, NULL);
    jdouble *Slm_imag = (*env)->GetDoubleArrayElements(env, imag, NULL);
    
    // Set the memory
    complex double *Slm = (complex double *) shtns_malloc( shtns->nlm * sizeof(complex double));
    double *Th = (double *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(double));
    
    // Set the values of Slm
    int index = 0;
    for (int m = 0; m <= shtns->mmax; m++) {
        for (int l = m; l <= shtns->lmax; l++) {
            Slm[index] = (double) Slm_real[index] + _Complex_I*((double) Slm_imag[index]);
            index++;
        }
    }
    
    // come back to spatial space: Slm -> Th
    SH_to_spat(shtns, Slm, Th);			// this does not destroy the spectral data in Slm.
    
    // Setup the return
    int arraySize = shtns->nphi * shtns->nlat;
    jdoubleArray spatialArray = (*env)->NewDoubleArray(env, arraySize);         // allocate
    (*env)->SetDoubleArrayRegion(env, spatialArray, 0 , arraySize, Th);   // copy
    
    // Cleanup
    (*env)->ReleaseDoubleArrayElements(env, real, Slm_real, 0);
    (*env)->ReleaseDoubleArrayElements(env, imag, Slm_imag, 0);
    return spatialArray;
}
