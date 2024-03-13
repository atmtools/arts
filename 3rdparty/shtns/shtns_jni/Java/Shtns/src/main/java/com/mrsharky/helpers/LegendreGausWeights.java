package com.mrsharky.helpers;

import static com.mrsharky.helpers.Utilities.arange;
import static com.mrsharky.helpers.Utilities.linspace;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Arrays;

/**
 *
 * @author Julien Pierret (julien.pierret@gmail.com) based on Matlab code by Greg von Winckel.
 */
public class LegendreGausWeights {
    private int _N;
    private int _N1;
    private int _N2;
    private double _a;
    private double _b;
    private double[] _x;
    private double[] _w;
    
    public double[] GetWeights() {
        return _w;
    }
    
    public double[] GetValues() {
        return _x; 
   }
    
    public LegendreGausWeights(int n, double a, double b) throws Exception {
        if (n < 0) {
            throw new Exception ("Input argument n must be greater than 1: " + n);
        }
        this._N = n-1;
        this._N1 = _N + 1;
        this._N2 = _N + 2;
        this._a = a;
        this._b = b;
        
        double[] xu = linspace(a,b,_N1);
        double[] zeroToN1 = arange(0,_N1,1);
        
        // Calculate inital guess       
        double[] left = DoubleArray.Multiply(DoubleArray.Add(DoubleArray.Multiply(zeroToN1, 2.0)
                            , (1.0))
                        ,(Math.PI)/(2.0*_N+2.0));
        left = DoubleArray.Cos(left);
        
        double[] right = DoubleArray.Multiply(xu, Math.PI*_N/_N2);
        right = DoubleArray.Multiply(DoubleArray.Sin(right),0.27/_N1);
        
        double[] y = DoubleArray.Add(left, right);

        // Legendre-Gauss CosineOfMatrixValuesandermonde Matrix
        double[][] L = new double[_N1][_N2];
        double[] Lp =  new double[_N2];
        
        double epsilonCompare = DoubleArray.Max(DoubleArray.Abs(DoubleArray.Add(y,2)));
        while (epsilonCompare > Math.ulp(1.0)) {
            L = DoubleArray.SetColumn(L, 0, 1.0);
            L = DoubleArray.SetColumn(L,1,y);
            
            double[] l_loop  = arange(2,_N1+1.0,1);
            
            for (double k : l_loop) {
                int k_ = (int) Math.round(k);
                          
                double[] loopLeft = DoubleArray.Multiply(DoubleArray.Multiply(y,DoubleArray.GetColumn(L, (k_-1))),(2.0*k-1.0));
                double[] loopRight = DoubleArray.Multiply(DoubleArray.GetColumn(L,(((int)Math.round(k))-2)), (k-1.0));
                
                double[] columnValueToSet = DoubleArray.Multiply(DoubleArray.Add(loopLeft, DoubleArray.Multiply(loopRight,-1.0)),(1.0/k));
                L = DoubleArray.SetColumn(L,k_,columnValueToSet);
            }
            
            // Derivative of LGVM
            double[] lgvm_numerator = DoubleArray.Multiply(DoubleArray.Add(DoubleArray.GetColumn(L,_N1-1), 
                    DoubleArray.Multiply(
                            DoubleArray.Multiply(y, DoubleArray.GetColumn(L,_N2-1)), -1.0))
                    , _N2);            
            double[] lgvm_denominator = DoubleArray.Power(DoubleArray.Add(DoubleArray.Multiply(DoubleArray.Power(y,2.0),-1.0),1.0), -1.0);                   
            Lp = DoubleArray.Multiply (lgvm_numerator,lgvm_denominator);
            
            double[] y0 = Arrays.copyOf(y, y.length);
            
            y = DoubleArray.Add(y0, DoubleArray.Multiply(DoubleArray.Multiply(DoubleArray.GetColumn(L,_N2-1),-1.0),
            DoubleArray.Power(Lp,-1.0)));
            
            epsilonCompare = DoubleArray.Max(DoubleArray.Abs(DoubleArray.Add(y, DoubleArray.Multiply(y0, -1.0))));
        }
        
        // Linear map from [-1,1] to [a,b]
        double[] x_left  = DoubleArray.Multiply(DoubleArray.Add(DoubleArray.Multiply(y, -1.0), 1.0), a);
        double[] x_right = DoubleArray.Multiply(DoubleArray.Add(y,1.0), b);
        
        _x = DoubleArray.Multiply(DoubleArray.Add(x_left, x_right), 1.0/2.0);
        
        // Compute the weights
        double[] w_left  = DoubleArray.Add(DoubleArray.Multiply(DoubleArray.Power(y, 2.0), -1.0), 1.0);
        double[] w_right = DoubleArray.Power(Lp, 2.0);
        double w_scalar = (b-a)*Math.pow(((double)_N2)/((double)_N1),2.0);
        
        _w = DoubleArray.Multiply(DoubleArray.Power(DoubleArray.Multiply(w_left, w_right), -1.0), w_scalar);
    } 
      
    public static void main(String args[]) {
        try {
            LegendreGausWeights lgwt = new LegendreGausWeights(100, -1.0, 1.0);
                        
            double[] x = lgwt.GetValues();
            double[] w = lgwt.GetWeights();
            
            DoubleArray.Print(x);
            DoubleArray.Print(w);
           
        } catch (Exception ex) {
            Logger.getLogger(LegendreGausWeights.class.getName()).log(Level.SEVERE, null, ex);
        }  
    }
}
