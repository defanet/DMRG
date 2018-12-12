package dmrg;

import Tools.Matrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Hamilton {
    RealMatrix H = null, H_ = null, _H = null;
    /*initial spin matrix, Ref:[1] Eqs(2)(4)*/
    RealMatrix Sz = new Array2DRowRealMatrix(new double[][]{{0.5, 0}, {0, -0.5}});
    RealMatrix Sp = new Array2DRowRealMatrix(new double[][]{{0, 1}, {0, 0}});
    RealMatrix Sn = Sp.transpose();
    RealMatrix SzSz, SpSn, I;
    RealMatrix Sz_, Sp_, Sn_, _Sz, _Sp, _Sn;
    
    public Hamilton() {
        /*Hamilton 2 spin problem, Ref:[1] Eqs(12)*/
        SzSz = Matrix.kron(Sz, Sz);
        SpSn = Matrix.kron(Sp, Sn).add(Matrix.kron(Sn, Sp));
        H = SzSz.add(SpSn.scalarMultiply(0.5));
    }

    /* i --> index of particle counted from left*/
    public RealMatrix swapFromLeft(int i) {
        
        if (i < 2) return null;
        if (i == 2) H_ = H;
        else { 
            /*Hamilton mulitiple spin, n-->number spin, Ref:[1] Eqs(17)(18)(19)*/
            H_ = Matrix.kron(H_, Matrix.Identity(2));
            I = Matrix.Identity((int) Math.pow(2, i - 2));
            Sz_ = Matrix.kron(I, Sz);
            Sp_ = Matrix.kron(I, Sp);
            Sn_ = Matrix.kron(I, Sn);
            SzSz = Matrix.kron(Sz_, Sz);
            SpSn = (Matrix.kron(Sp_, Sn)).add(Matrix.kron(Sn_, Sp));
            H_ = H_.add(SzSz.add(SpSn.scalarMultiply(0.5)));
        }
        return H_;
    }
    
    public RealMatrix swapFromRight(int i) {
        
        if (i < 2) return null;
        if (i == 2) _H = H;
        else { 
            /*Hamilton mulitiple spin, n-->number spin, Ref:[1] Eqs(17)(18)(19)*/
            _H = Matrix.kron(Matrix.Identity(2), _H);
            I = Matrix.Identity((int) Math.pow(2, i - 2));
            _Sz = Matrix.kron(Sz, I);
            _Sp = Matrix.kron(Sp, I);
            _Sn = Matrix.kron(Sn, I);
            SzSz = Matrix.kron(Sz, _Sz);
            SpSn = (Matrix.kron(Sn, _Sp)).add(Matrix.kron(Sp, _Sn));
            _H = _H.add(SzSz.add(SpSn.scalarMultiply(0.5)));
        }
        return _H;
    }
    
    public RealMatrix combineSwap(RealMatrix H_, RealMatrix _H) {
        int dimL = H_.getRowDimension();
        int dimR = _H.getRowDimension();
        H = Matrix.kron(H_, Matrix.Identity(dimR)).add(Matrix.kron(Matrix.Identity(dimR),_H));
        I = Matrix.Identity(dimL/2);
        Sz_ = Matrix.kron(I, Sz);
        Sp_ = Matrix.kron(I, Sp);
        Sn_ = Matrix.kron(I, Sn);
        _Sz = Matrix.kron(Sz, I);
        _Sp = Matrix.kron(Sp, I);
        _Sn = Matrix.kron(Sn, I);
        System.out.println(dimL+"\t"+dimR+"\t"+H.getRowDimension()+"\t"+Sz_.getRowDimension());
        H = H.add(Matrix.kron(Sz_, _Sz));
        H = H.add(Matrix.kron(Sp_, _Sn).scalarMultiply(0.5));
        H = H.add(Matrix.kron(Sn_, _Sp).scalarMultiply(0.5));
        return H;
    }
    
    public void print_blocks(int l,int r) {
        System.out.println("**********************************");
        System.out.println("LEFT SIZE = "+l);
        System.out.println("RIGHT SIZE = "+r);
        while (--l > 0) System.out.print("X ");
        if(r>0) System.out.print("* * ");
        else System.out.print("* ");
        while (--r > 0) System.out.print("X ");
        System.out.println();
    }
}
