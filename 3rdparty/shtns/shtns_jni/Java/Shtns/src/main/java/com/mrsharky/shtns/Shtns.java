package com.mrsharky.shtns;

import com.mrsharky.helpers.LegendreGausWeights;
import org.apache.commons.math3.complex.Complex;
import com.mrsharky.helpers.DoubleArray;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Julien Pierret (julien.pierret@gmail.com)
 */
public class Shtns {
    
    // Native Methods
    private native long initializeGaussian(int Q, int Ql, int periodicity, int jnlat, int jnlong, int numThreads);
    private native double[] spatialToSpectral(long cfg, double[] data);
    private native double[] spectralToSpatial(long cfg, double[] real, double[] imag);
    private native long getNlm(long cfg);

    /**
    * Adds the specified path to the java library path
    * From: http://fahdshariff.blogspot.be/2011/08/changing-java-library-path-at-runtime.html
    * @param pathToAdd the path to add
    * @throws Exception
    */
    public static void addLibraryPath(String pathToAdd) throws Exception{
        final Field usrPathsField = ClassLoader.class.getDeclaredField("usr_paths");
        usrPathsField.setAccessible(true);

        //get array of paths
        final String[] paths = (String[])usrPathsField.get(null);

        //check if the path to add is already present
        for(String path : paths) {
            if(path.equals(pathToAdd)) {
                return;
            }
        }

        //add the new path
        final String[] newPaths = Arrays.copyOf(paths, paths.length + 1);
        newPaths[newPaths.length-1] = pathToAdd;
        usrPathsField.set(null, newPaths);
    }
    
    // Load the necessary native library
    static {
        try {
            addLibraryPath("/usr/local/lib");
            addLibraryPath("/media/dropbox/PhD/Reboot/Projects/Ncep20thCenturyReanalysisV2c/C/shtns_jni/C/dist/Debug/GNU-Linux");
            System.loadLibrary("shtns_jni"); // Load native library at runtime
            //System.load("/media/dropbox/PhD/Reboot/Projects/Ncep20thCenturyReanalysisV2c/C/shtns_jni/dist/Debug/GNU-Linux/libshtns_jni.so");                                       
        } catch (Exception ex) {
            Logger.getLogger(Shtns.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private final long _cfg;
    private final int _lmax;
    private final int _mmax;
    private final int _mres;
    private final int _nlat;
    private final int _nphi;
    private final long _nlm;
    private Complex[] _spectral;
    private double[] _spatial;
    
    public Shtns(int lmax, int mmax, int mres, int nlat, int nphi, int numThreads) throws Exception {
        _lmax = lmax;
        _mmax = mmax;
        _mres = mres;
        _nlat = nlat;
        _nphi = nphi;
        _cfg = initializeGaussian(_lmax, _mmax, _mres, _nlat, _nphi, numThreads);
        _nlm = getNlm(_cfg);
    }
    
    /**
     * 
     * @param data Complex Array expected to be in the "compressed" format
     * @return 
     */
    public void SpectralToSpatial(Complex[] data) {
        _spectral = data;
        // Split the real and imag apart
        double[] real = new double[data.length];
        double[] imag = new double[data.length];
        
        for (int i = 0; i < data.length; i++) {
            real[i] = data[i].getReal();
            imag[i] = data[i].getImaginary();
        }   
        
        double[] jniData = spectralToSpatial(_cfg, real, imag);
        _spatial = jniData;
    }
    
    public double[][] GetSpatial() {
        double[][] results = new double[_nlat][_nphi];
        for (int ip = 0; ip < _nphi; ip++) {		// loop on longitude
            for (int it = 0; it < _nlat; it++) {	// loop on latitude
                results[it][ip] = _spatial[ip*_nlat + it];
            }
	}
        return results;
    }
    
    public Complex[] GetSpectralCompressed() throws Exception {
        if (this._spectral != null) {
            return this._spectral;
        } else {
            throw new Exception("Spectral has not been processed!");
        }
    }
    
    public double[] GetSpatialCompressed() throws Exception {
        if (this._spatial != null) {
            return this._spatial;
        } else {
            throw new Exception("Spatial has not been processed!");
        }
    }
    
    /**
     * 
     * @param data Array must be organized by long (- to +) & lat (- to +)
     */
    public void SpatialToSpectral(double[] data) throws Exception {
        _spatial = data;
        double[] jniData = spatialToSpectral(_cfg, data);
        if (jniData.length != (this._nlm*2)) {
            throw new Exception("Returned spectral length from SHtns doesn't match expected length");
        }
        Complex[] results = new Complex[(int) this._nlm];
    
        int counter = 0;
        for (int m = 0; m <= _mmax; m++) {
            for (int l = m; l <= _lmax; l++) {
                int index = counter;
                double real = jniData[2*index];
                double imag = jniData[2*index+1];
                results[index] = new Complex(real, imag);
                counter++;
            }
        }
        this._spectral = results;
    }
    
    /**
     * 
     * @param data Double array must be organized double[south to north][east to west]
     * @return
     * @throws Exception 
     */
    public void SpatialToSpectral(double[][] data) throws Exception {
        
        // double check the dimensions of the new data are correct
        if (data.length != _nlat) {
            throw new Exception("Input data number of latitude points don't match what was configured");
        }
        if (data[0].length != _nphi) {
            throw new Exception("Input data number of longitude points don't match what was configured");
        }
        
        // first re-organize the data into a single double array
        double[] dataToShtns = new double[_nphi*_nlat];
        
        for (int ip = 0; ip < _nphi; ip++) {		// loop on longitude
            for (int it = 0; it < _nlat; it++) {	// loop on latitude
                dataToShtns[ip*_nlat + it] = data[it][ip];
            }
        }
        SpatialToSpectral(dataToShtns);
    }
    
    // Test Driver
    public static void main(String[] args) throws Exception { 

        int lmax = 5;		// maximum degree of spherical harmonics
        int mmax = 3;		// maximum order of spherical harmonics
        int mres = 1;		// periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)
        int nlat = 32;		// number of points in the latitude direction  (constraint: nlat >= lmax+1)
        int nphi = 10;        

        LegendreGausWeights lgw = new LegendreGausWeights(nlat,-1,1);
        
        double[] legZerosM = lgw.GetValues();
        double[] legZerosRad = DoubleArray.Add(DoubleArray.ArcSin(legZerosM), (Math.PI/2.0));
        double[] phi2 = DoubleArray.Cos(legZerosRad);

        double[][] data = new double[nlat][nphi];
        double[] data2 = new double[nlat*nphi];

        // generate some gridded data:
        for (int ip = 0; ip < nphi; ip++) {		// loop on longitude
            double phi = (2.*Math.PI/(mres*nphi)) * ip;
            double cos_phi = Math.cos(phi);
            for (int it = 0; it < nlat; it++) {		// loop on latitude
                double cos_theta = phi2[it];
                double sin_theta = Math.sqrt(1.0 - cos_theta*cos_theta);
                double currValue = cos_phi*cos_theta*sin_theta;     // this is an m=1, l=2 harmonic
                data[it][ip] = currValue;
                data2[ip*nlat + it] = currValue;
            }
        }
        
        // Setup shtns
        System.out.println("Setting up SHTns");
        Shtns shtns = new Shtns(lmax, mmax, mres, nlat, nphi, 0);
        
        // spatial to spectral
        System.out.println("Spatial to Spectral");
        shtns.SpatialToSpectral(data2);
        Complex [] spectral = shtns.GetSpectralCompressed();
        
        // spectral back to spatial
        System.out.println("Spectral to Spatial");
        shtns.SpectralToSpatial(spectral);
        double[][] dataRebuilt = shtns.GetSpatial();
        
        DoubleArray.Print(data);
        DoubleArray.Print(dataRebuilt);
        
        // Compare the two
        double error = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                error = error + Math.pow(data[i][j] - dataRebuilt[i][j], 2.0);
            }
        }
        System.out.println("SumSqError: " + error);  
    }
}