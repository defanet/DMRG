package dmrg;

import Tools.Matrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * @author defa, Nov 14, 2018
 * Reference:
 * [1] Adrian E. Feiguin, The Density Matrix Renormalization Group and its time-dependent variants
 */
public class HeisenbergChain {
    GroundState groundState = new GroundState();
    /*initial spin matrix, Ref:[1] Eqs(2)(4)*/
    RealMatrix Sz = new Array2DRowRealMatrix(new double[][] {{0.5, 0},{0, -0.5}});
    RealMatrix Sp = new Array2DRowRealMatrix(new double[][] {{0, 1},{0, 0}});
    RealMatrix Sn = Sp.transpose();
    
    /*Hamilton mulitiple spin, n-->number spin, Ref:[1] Eqs(17)(18)(19)*/
    public RealMatrix Hamilton(int n) {
        RealMatrix H = null;
        RealMatrix Sz_, Sp_, Sn_, SzSz, SpSn, I;
        if(n<2) return null; 

        for(int i=2; i<=n; i++) {
            if(i==2) {
                /*Hamilton 2 spin problem, Ref:[1] Eqs(12)*/
                SzSz = Matrix.kron(Sz, Sz);
                SpSn = Matrix.kron(Sp, Sn).add(Matrix.kron(Sn, Sp));
                H = SzSz.add(SpSn.scalarMultiply(0.5));
            }
            else {
                H = Matrix.kron(H, MatrixUtils.createRealIdentityMatrix(2));
                I = MatrixUtils.createRealIdentityMatrix((int)Math.pow(2,i-2));
                Sz_ = Matrix.kron(I, Sz);
                Sp_ = Matrix.kron(I, Sp);
                Sn_ = Matrix.kron(I, Sn);            
                SzSz = Matrix.kron(Sz_, Sz);
                SpSn = (Matrix.kron(Sp_, Sn)).add(Matrix.kron(Sn_, Sp));
                H = H.add(SzSz.add(SpSn.scalarMultiply(0.5)));
            }
            groundState.lanczos(H, null, true);
        }
        return H;
    }
    
    public static void main(String args[]) {
        int numParticle = 6;
        HeisenbergChain Hs = new HeisenbergChain();
        RealMatrix m = Hs.Hamilton(numParticle);
        /*double[][] triDiagMatrix = getTridiagonalMatrix(m.getData(), 23);
        Matrix.print(triDiagMatrix);
        EigenDecomposition eigen = new EigenDecomposition(new Array2DRowRealMatrix(triDiagMatrix));
        double[] eigenValue = eigen.getRealEigenvalues();
        Array.print(eigenValue);
        System.out.println("--------------------");
        Array.print(Array.shortDesc(eigenValue));*/
    }
}
