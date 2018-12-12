package Tools;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Matrix extends Array2DRowRealMatrix {
    
    public Matrix(double[][] m) {super(m);}
    public Matrix(int m, int n) {super(m,n);}
    
    public static RealMatrix kron(RealMatrix a, RealMatrix b) {
        int aRow = a.getRowDimension();
        int aCol = a.getColumnDimension();
        int bRow = b.getRowDimension();
        int bCol = b.getColumnDimension();
        if(aRow==1) return b;
        if(bRow==1) return a;
        RealMatrix m = new Array2DRowRealMatrix(aRow*bRow, aCol*bCol);
        for(int i=0; i<aRow; i++) 
            for(int k=0; k<bRow; k++) 
                for(int j=0; j<aCol; j++) 
                    for(int l=0; l<bCol; l++) 
                        m.setEntry(i*bRow+k, j*bCol+l, a.getEntry(i, j)*b.getEntry(k, l));
        return m;
    }
    
    public static RealMatrix Identity(int dim) {
        return MatrixUtils.createRealIdentityMatrix(dim);
    }
    
    public static RealMatrix random(int m, int n) {
        RealMatrix matrix = new Array2DRowRealMatrix(m, n);
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                matrix.setEntry(i, j, Math.random());
        return matrix;
    }
    
    
    public static void print(RealMatrix m) {
        System.out.println();
        for(int i=0; i<m.getRowDimension(); i++) {
            //System.out.print("[\t");
            for(int j=0; j<m.getColumnDimension(); j++) 
                System.out.printf("%8.2f ",m.getEntry(i, j));
            System.out.println();
        }
    }
    
    public static void print(String str, RealMatrix m) {System.out.println(str);print(m);}
    
    public static void print(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%8.2f ", matrix[i][j]);
            }
            System.out.println();
        }
    }
    
    public static void print(RealVector v) {
        for (int i = 0; i < v.getDimension(); i++)
                System.out.print(v.getEntry(i)+"\t");
        System.out.println();
    }
    
    public static void print(String str, RealVector v) {System.out.println(str); print(v);}
    
}
