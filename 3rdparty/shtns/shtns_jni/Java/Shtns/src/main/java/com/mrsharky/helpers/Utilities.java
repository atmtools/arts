package com.mrsharky.helpers;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.ArrayFieldVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.FieldVector;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * @author Julien Pierret (julien.pierret@gmail.com)
 */
public class Utilities {
    
    public static double[] ArrayOnes(int size) {
        double[] ret = new double[size];
        for (int counter = 0; counter < ret.length; counter++) {
            ret[counter] = 1;
        } 
        return ret;
    }
    
    
    public static int[] minToMaxByOne(int min, int max) {
        int[] ret = new int[max-min+1];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = min + i;
        }
        return ret;
    }
    
    public static double[] minToMaxByOne(double min, double max) {
        double[] ret = new double[(int) max - (int) min+1];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = min + i;
        }
        return ret;
    }
    
    public static double[] linspace(double min, double max, int points) {  
        double[] d = new double[points]; 
        if (points == 1) {
            d[0] = min;
        } else {
            for (int i = 0; i < points; i++){  
                d[i] = min + i * (max - min) / (points - 1);  
            }
        }
        return d;  
    }
    
    public static int[] linspace(int min, int max, int points) {  
        int[] d = new int[points]; 
        if (points == 1) {
            d[0] = min;
        } else {
            for (int i = 0; i < points; i++){  
                d[i] = min + i * (max - min) / (points - 1);  
            }
        }
        return d;  
    }
    
    public static double[] arange(double start, double end, double step) {
        int points = (int) Math.floor((end-start)/step);
        double[] d = new double[points];
        for (int i = 0; i < points; i++) {
            d[i] = start + step*i;
        }
        return d;
    }
    
    public static RealVector CosineOfVectorValues (RealVector vec) {
        RealVector ret = vec.copy();
        for (int counter = 0; counter < ret.getDimension(); counter++) {
            ret.setEntry(counter, Math.cos(ret.getEntry(counter)));
        }
        return ret;
    }
    
    public static RealVector SineOfVectorValues (RealVector vec) {
        RealVector ret = vec.copy();
        for (int counter = 0; counter < ret.getDimension(); counter++) {
            ret.setEntry(counter, Math.sin(ret.getEntry(counter)));
        }
        return ret;
    }
    
    
    public static int ComplexExample() {
        // create a 2x2 complex matrix 
        Complex[][] matrixData = new Complex[][] { 
            { new Complex(1.0,  0.0), new Complex( 0.0, 1.0) }, 
            { new Complex(0.0, -1.0), new Complex(-1.0, 0.0) } 
        }; 
        FieldMatrix<Complex> m = new Array2DRowFieldMatrix<Complex>(matrixData); 


        // create a vector 
        Complex[] vectorData = new Complex[] { 
            new Complex(1.0, 2.0), 
            new Complex(3.0, 4.0), 
        }; 
        FieldVector<Complex> u = new ArrayFieldVector<Complex>(vectorData); 

        // perform matrix-vector multiplication 
        FieldVector<Complex> v = m.operate(u); 
        

        // print the initial vector 
        for (int i = 0; i < u.getDimension(); ++i) { 
            //System.out.println(ComplexFormat.formatComplex(u.getEntry(i))); 
        } 

        System.out.println(); 

        // print the result 
        for (int i = 0; i < v.getDimension(); ++i) { 
            //System.out.println(ComplexFormat.formatComplex(v.getEntry(i))); 
        } 

        return 1;
    }
}
