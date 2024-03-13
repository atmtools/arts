package com.mrsharky.helpers;

/**
 *
 * @author Julien Pierret (julien.pierret@gmail.com)
 */
public class DoubleArray {

    public static double[] Multiply(double[] mat, double scaler) {
        double[] ret = new double[mat.length];
        for (int row = 0; row < mat.length; row++) {
            ret[row] = mat[row] * scaler;
        }
        return ret;
    }

    public static double[] Multiply(double[] mat1, double[] mat2) throws Exception {
        if (mat1.length != mat2.length) {
            throw new Exception("Array dimensions do not agree");
        }
        double[] ret = new double[mat1.length];
        for (int i = 0; i < mat1.length; i++) {
            ret[i] = mat1[i] * mat2[i];
        }
        return ret;
    }

    public static double[] Cos(double[] mat) {
        double[] ret = new double[mat.length];
        for (int i = 0; i < mat.length; i++) {
            ret[i] = Math.cos(mat[i]);
        }
        return ret;
    }

    public static double Max(double[] mat) {
        double maxValue = -Double.MAX_VALUE;
        for (int row = 0; row < mat.length; row++) {
            double currDouble = mat[row];
            if (currDouble > maxValue) {
                maxValue = currDouble;
            }
        }
        return maxValue;
    }

    public static double[][] SetColumn(double[][] mat, int column, double value) {
        double[][] ret = new double[mat.length][mat[0].length];
        for (int row = 0; row < mat.length; row++) {
            for (int col = 0; col < mat[row].length; col++) {
                if (col == column) {
                    ret[row][col] = value;
                } else {
                    ret[row][col] = mat[row][col];
                }
            }
        }
        return ret;
    }

    public static double[][] SetColumn(double[][] mat, int column, double[] value) {
        double[][] ret = new double[mat.length][mat[0].length];
        for (int row = 0; row < mat.length; row++) {
            for (int col = 0; col < mat[row].length; col++) {
                if (col == column) {
                    ret[row][col] = value[row];
                } else {
                    ret[row][col] = mat[row][col];
                }
            }
        }
        return ret;
    }
    
    public static double[] ArcSin (double[] mat) {
        double[] output = new double[mat.length];
        for (int i = 0; i < mat.length; i++) {
            output[i] = Math.asin(mat[i]);
        }
        return output;
    }

    public static void Print(double[][] values) {
        System.out.println();
        for (int col = 0; col < values[0].length; col++) {
            System.out.print("\t" + (col + 1));
        }
        System.out.println();
        for (int row = 0; row < values.length; row++) {
            System.out.print((row + 1) + ":\t");
            for (int col = 0; col < values[0].length; col++) {
                System.out.print(values[row][col] + "\t");
            }
            System.out.println();
        }
    }

    public static void Print(double[] values) {
        System.out.println();
        System.out.println("\t" + 1);
        for (int row = 0; row < values.length; row++) {
            System.out.print((row + 1) + ":\t");
            System.out.println(values[row] + "\t");
        }
    }

    public static double[] Sin(double[] mat) {
        double[] ret = new double[mat.length];
        for (int row = 0; row < mat.length; row++) {
            ret[row] = Math.sin(mat[row]);
        }
        return ret;
    }

    public static double[] Power(double[] mat, double exponent) {
        double[] ret = new double[mat.length];
        for (int row = 0; row < mat.length; row++) {
            double val = mat[row];
            ret[row] = Math.pow(val, exponent);
        }
        return ret;
    }

    public static double[] GetColumn(double[][] mat1, int column) {
        double[] ret = new double[mat1.length];
        for (int row = 0; row < mat1.length; row++) {
            ret[row] = mat1[row][column];
        }
        return ret;
    }

    public static double[] Add(double[] mat, double scaler) {
        double[] ret = new double[mat.length];
        for (int row = 0; row < mat.length; row++) {
            ret[row] = mat[row] + scaler;
        }
        return ret;
    }
    
    public static double[] Add(double[] mat1, double[] mat2) throws Exception {
        if (mat1.length != mat2.length) {
            throw new Exception("Array dimensions do not agree");
        }
        double[] ret = new double[mat1.length];
        for (int i = 0; i < mat1.length; i++) {
            ret[i] = mat1[i] + mat2[i];
        }
        return ret;
    }

    public static double[] Abs(double[] mat) {
        double[] ret = new double[mat.length];
        for (int row = 0; row < mat.length; row++) {
            ret[row] = Math.abs(mat[row]);
        }
        return ret;
    }    
}
