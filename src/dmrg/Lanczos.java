package dmrg;

import Tools.Matrix;
import Tools.Array;
import org.apache.commons.math3.linear.*;

public class Lanczos {
    
    

    public static double[][] getTridiagonalMatrix(double[][] a, int dim) {
        double[] alpha = new double[dim + 1];
        double[] betta = new double[dim + 1];

        int n = a.length;

        RealVector v0 = new ArrayRealVector(n);
        RealVector v1 = new ArrayRealVector(Array.getRandom(n));

        RealVector wx = new ArrayRealVector(n);
        RealVector w = new ArrayRealVector(n);

        RealMatrix t = new Array2DRowRealMatrix(a);

        for (int i = 1; i < dim; i++) {
            wx = t.operate(v1);
            alpha[i] = wx.dotProduct(v1);
            w = wx.subtract(v1.mapMultiply(alpha[i])).subtract(v0.mapMultiply(betta[i]));
            betta[i + 1] = w.getNorm();
            v0 = v1.copy();
            v1 = w.mapMultiply(1 / betta[i + 1]);
        }

        return makeTridiagonalMatrix(alpha, betta);
    }

    public static double[][] makeTridiagonalMatrix(double[] alpha, double[] betta) {
        int m = alpha.length - 1;
        double[][] T = new double[m][m];

        T[0][0] = alpha[1];
        for (int i = 1; i < alpha.length - 1; i++) {
            T[i][i] = alpha[i + 1];
            T[i - 1][i] = betta[i + 1];
            T[i][i - 1] = betta[i + 1];
        }

        return T;
    }

    

    public static void main(String[] args) {
        HeisenbergChain Hs = new HeisenbergChain();
        RealMatrix m = Hs.Hamilton(10);
        double[][] triDiagMatrix = getTridiagonalMatrix(m.getData(), 23);
        Matrix.print(triDiagMatrix);
        EigenDecomposition eigen = new EigenDecomposition(new Array2DRowRealMatrix(triDiagMatrix));
        double[] eigenValue = eigen.getRealEigenvalues();
        Array.print(eigenValue);
        System.out.println("--------------------");
        Array.print(Array.shortDesc(eigenValue));
    }
}
