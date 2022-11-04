package com.mrsharky.helpers;

import org.apache.commons.math3.complex.Complex;

/**
 *
 * @author Julien Pierret (julien.pierret@gmail.com)
 */
public class ComplexArray {
   
    public static Complex[] Log (Complex[] comp1) throws Exception {
        Complex[] comp = new Complex[comp1.length];
        for (int i = 0; i < comp1.length; i++) {
            comp[i] = comp1[i].log();
        }
        return comp;
    }
    
    public static Complex[][] Log (Complex[][] comp1) throws Exception {
        Complex[][] comp = new Complex[comp1.length][comp1[0].length];
        for (int i = 0; i < comp1.length; i++) {
            for (int j = 0; j < comp1[0].length; j++) {
                comp[i][j] = comp1[i][j].log();
            }
        }
        return comp;
    }
    
    public static Complex[] Log (double[] comp1) throws Exception {
        Complex[] comp = new Complex[comp1.length];
        for (int i = 0; i < comp1.length; i++) {
            Complex currValue = new Complex(comp1[i]);
            comp[i] = currValue.log();
        }
        return comp;
    }
    
    public static Complex[][] Log (double[][] comp1) throws Exception {
        Complex[][] comp = new Complex[comp1.length][comp1[0].length];
        for (int i = 0; i < comp1.length; i++) {
            for (int j = 0; j < comp1[0].length; j++) {
                Complex currValue = new Complex(comp1[i][j]);
                comp[i][j] = currValue.log();
            }
        }
        return comp;
    }
    
    public static Complex[] Multiply (Complex[] comp1, double comp2) throws Exception {
        Complex[] comp = new Complex[comp1.length];
        for (int i = 0; i < comp1.length; i++) {
            comp[i] = comp1[i].multiply(comp2);
        }
        return comp;
    }
    
    public static Complex[][] Multiply (Complex[][] comp1, double comp2) throws Exception {
        Complex[][] comp = new Complex[comp1.length][comp1[0].length];
        for (int i = 0; i < comp1.length; i++) {
            for (int j = 0; j < comp1[0].length; j++) {
                comp[i][j] = comp1[i][j].multiply(comp2);
            }
        }
        return comp;
    }
    
    public static Complex[][] Multiply (Complex[][] comp1, Complex comp2) throws Exception {
        Complex[][] comp = new Complex[comp1.length][comp1[0].length];
        for (int i = 0; i < comp1.length; i++) {
            for (int j = 0; j < comp1[0].length; j++) {
                comp[i][j] = comp1[i][j].multiply(comp2);
            }
        }
        return comp;
    }
    
    public static Complex[][] Conjugate (Complex[][] comp1) throws Exception {
        Complex[][] comp = new Complex[comp1.length][comp1[0].length];
        for (int i = 0; i < comp1.length; i++) {
            for (int j = 0; j < comp1[0].length; j++) {
                comp[i][j] = comp1[i][j].conjugate();
            }
        }
        return comp;
    }
    
    public static double[][] Real (Complex[][] comp1) throws Exception {
        double[][] comp = new double[comp1.length][comp1[0].length];
        for (int i = 0; i < comp1.length; i++) {
            for (int j = 0; j < comp1[0].length; j++) {
                comp[i][j] = comp1[i][j].getReal();
            }
        }
        return comp;
    }
    
    public static void Print(Complex[][][] values) {
        System.out.println();
        for (int set = 0; set < values[0][0].length; set++) {
            System.out.println("Set=" + set);
            for (int col = 0; col <values[0].length; col++) { 
                System.out.print("\t" + (col+1));
            }
            System.out.println();
            for (int row = 0; row < values.length; row++) {
                System.out.print((row+1) +":\t");
                for (int col = 0; col <values[0].length; col++) {
                    System.out.print(values[row][col][set] + "\t");
                }
                System.out.println();
            }
        }
    }
    
    public static void Print(Complex[][] values) {
        System.out.println();
        for (int col = 0; col < values[0].length; col++) {
            System.out.print("\t" + (col + 1));
        }
        System.out.println();
        for (int row = 0; row < values.length; row++) {
            System.out.print((row + 1) + ":\t");
            for (int col = 0; col < values[0].length; col++) {
                Complex curr = values[row][col];
                System.out.print("(" + curr.getReal() + ", " + curr.getImaginary() +  ")\t");
            }
            System.out.println();
        }
    }
    
    public static void Print(Complex[] values) {
        System.out.println();
        for (int col = 0; col < 1; col++) {
            System.out.print("\t" + (col + 1));
        }
        System.out.println();
        for (int row = 0; row < values.length; row++) {
            System.out.print((row + 1) + ":\t");
            Complex curr = values[row];
            System.out.print("(" + curr.getReal() + ", " + curr.getImaginary() +  ")\n");
        }
        System.out.println();
    }
    
    public static Complex[] CreateComplex (int vec1) {
        Complex[] comp = new Complex[vec1];
        for (int i = 0; i < vec1; i++) {
            comp[i] = new Complex(0.0,0.0);
        }
        return comp;
    }
    
    public static Complex[] CreateComplex (double[] vec1) {
        Complex[] comp = new Complex[vec1.length];
        for (int i = 0; i < vec1.length; i++) {
            comp[i] = new Complex(vec1[i],0.0);
        }
        return comp;
    }
    
    public static Complex[][] CreateComplex (int vec1, int vec2) {
        Complex[][] comp = new Complex[vec1][vec2];
        for (int i = 0; i < vec1; i++) {
            for (int j = 0; j < vec2; j++) {
                comp[i][j] = new Complex(0.0,0.0);
            }
        }
        return comp;
    }
    
    public static Complex[][][] CreateComplex (int vec1, int vec2, int vec3) {
        Complex[][][] comp = new Complex[vec1][vec2][vec3];
        for (int i = 0; i < vec1; i++) {
            for (int j = 0; j < vec2; j++) {
                for (int k = 0; k < vec3; k++) {
                comp[i][j][k] = new Complex(0.0,0.0);
                }
            }
        }
        return comp;
    }
    
    public static Complex[] RowSum(Complex[][] origData) {
        Complex[] ret = new Complex[origData.length];
        for (int row = 0; row < origData.length; row++) {
            Complex total = new Complex(0.0, 0.0);
            for (int col = 0; col < origData[0].length; col++) {
                total = total.add(origData[row][col]);
            }
            ret[row] = total;
        }
        return ret;
    }
    
    public static Complex Sum(Complex[] origData) {
        Complex total = new Complex(0,0);
        for (int row = 0; row < origData.length; row++) {
            total = total.add(origData[row]);
        }
        return total;
    }
    
    public static Complex[][][] SetData(Complex[][][] origData, Complex[][] newData, int[] vec1, int[] vec2, int vec3) {
        for (int vec1Counter = 0; vec1Counter < vec1.length; vec1Counter++) {
            int currVec1 = vec1[vec1Counter];
            for (int vec2Counter = 0; vec2Counter < vec2.length; vec2Counter++) {
                int currVec2 = vec2[vec2Counter];
                origData[currVec1][currVec2][vec3] = newData[vec1Counter][vec2Counter];
            }
        }
        return origData;
    }
    
    public static Complex[][] SetData(Complex[][] origData, Complex[] newData, int vec1, int[] vec2) {
        for (int vec2Counter = 0; vec2Counter < vec2.length; vec2Counter++) {
            int currVec2 = vec2[vec2Counter];
            origData[vec1][currVec2] = newData[vec2Counter];
        }
        return origData;
    }
    
    public static Complex[][] GetData(Complex[][][] origData, int[] vec1, int[] vec2, int vec3) {
        Complex[][] newData = new Complex[vec1.length][vec2.length];
        for (int vec1Counter = 0; vec1Counter < vec1.length; vec1Counter++) {
            int currVec1 = vec1[vec1Counter];
            for (int vec2Counter = 0; vec2Counter < vec2.length; vec2Counter++) {
                int currVec2 = vec2[vec2Counter];
                newData[vec1Counter][vec2Counter] = origData[currVec1][currVec2][vec3];
            }
        }
        return newData;
    }
    
    public static Complex[][] GetData(Complex[][][] origData, int vec1, int[] vec2, int[] vec3) {
        Complex[][] newData = new Complex[vec2.length][vec3.length];
        for (int vec2Counter = 0; vec2Counter < vec2.length; vec2Counter++) {
            int currVec2 = vec2[vec2Counter];
            for (int vec3Counter = 0; vec3Counter < vec3.length; vec3Counter++) {
                int currVec3 = vec3[vec3Counter];
                newData[vec2Counter][vec3Counter] = origData[vec1][currVec2][currVec3];
            }
        }
        return newData;
    }
    
    public static Complex[] GetData(Complex[][] origData, int vec1, int[] vec2) {
        Complex[] newData = new Complex[vec2.length];
        for (int vec2Counter = 0; vec2Counter < vec2.length; vec2Counter++) {
            int currVec2 = vec2[vec2Counter];
            newData[vec2Counter] = origData[vec1][currVec2];
        }
        return newData;
    }
   
}
