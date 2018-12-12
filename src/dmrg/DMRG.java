package dmrg;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * @author defa, Nov 14, 2018
 * Reference:
 * [1] Adrian E. Feiguin, The Density Matrix Renormalization Group and its time-dependent variants
 */
public class DMRG {
    RealMatrix H, H_, _H;
    Hamilton block;
    GroundState groundState;
    
    public DMRG(int n) {
        int dim = 0;
        block = new Hamilton();
        groundState = new GroundState();
        for(int i=2; i<=n/2; i++) {
            block.print_blocks(i, i);
            H_ = block.swapFromLeft(i);
            _H = block.swapFromRight(i);
        }
        H = block.combineSwap(H_, _H);
        //Matrix.print(H);
        groundState.lanczos(H, null, true);
    }
    
    
    public static void main(String args[]) {
        int numParticle = 12;
        DMRG dmrg = new DMRG(numParticle);
    }
}
