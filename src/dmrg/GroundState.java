package dmrg;

import Tools.Array;
import java.util.ArrayList;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * @author defa Copyright Lanzcos subrutine C++
 */
public class GroundState {

    private double SIGN(double a, double b) {
        if (b < 0) {
            return -Math.abs(a);
        } else {
            return Math.abs(a);
        }
    }

    /*
void lanczos(1. RealMatrix H,
             2. RealVector gs, 3. double[] seed, 
             4. double ener,
             5. double[] a, 6. double[] b, 
             7. int maxiter, 8. double tol,
             9. boolean use_seed, 10. boolean calc_gs, 11. String vectors_file, 12. boolean force_maxiter)    
     */
    public void lanczos(RealMatrix H, double[] seed, boolean calc_gs) {
        int dim = H.getRowDimension();
        double[] a = new double[1000];
        double[] b = new double[1000];
        boolean force_maxiter = false;
        int maxiter = -1;
        double ener;
        double tolerance = 1.e-7;
        ArrayList<RealVector> save = new ArrayList<>();
        RealVector q;
        RealVector gs = new ArrayRealVector(dim);
        RealVector x2 = new ArrayRealVector(dim);
        RealVector aux = new ArrayRealVector(dim);

        double[] d = new double[1000];
        double[] e = new double[1000];
        RealVector xc;
        double[][] z = new double[1000][1000];
        int col = 0, iter;
        double eini, e0;
        int control_max = maxiter;

        if (maxiter == -1) {
            force_maxiter = false;
        }

        //if(control_max == 0) { gs = T.get(1); maxiter = 1; return; }
        e0 = 10000.0f;
        maxiter = 0;
        //the fact, this while iteration only run one time
        while (true) // Two iterations: 1) calculate energy 2) build gs
        {
            Array.setValues(a, 0.0f);
            Array.setValues(b, 0.0f);

            if (seed != null) {
                q = new ArrayRealVector(seed);
            } else {
                q = new ArrayRealVector(Array.getRandom(dim));
            }

            b[0] = q.getNorm();
            q = q.mapDivide(b[0]);
            x2.set(0.0);
            b[0] = 1.;

            int nmax = Math.min(999, dim);
            for (iter = 1; iter <= nmax; iter++) { // Lanczos iteration
                eini = e0;

                if (b[iter - 1] != 0.) {
                    aux = q;
                    q = x2.mapMultiply(-b[iter - 1]);
                    x2 = aux.mapDivide(b[iter - 1]);
                    //Array.print(q.toArray());
                }
                aux.set(0.0);
                aux = H.operate(x2);
                q = q.add(aux);
                a[iter] = aux.dotProduct(x2);
                q = q.subtract(x2.mapMultiply(a[iter]));
                b[iter] = q.getNorm();
                //System.out.println(b[iter]);

                save.add(x2);

                //cout << "Iter=" << iter << " " << a(iter) << " " << b(iter) << endl;
                if (maxiter > 0) {  // We are building the ground state;
                } else {  // we are calculating energy

                    if (iter >= 2) {
                        Array.copy(d, 1, iter, a, 1, iter);
                        Array.copy(e, 2, iter + 1, b, 1, iter);
                        // call tqli without eigenvectors
                        tqli(d, e, iter, z, false);

                        e0 = 10000.0f;
                        for (int j = 1; j <= iter; j++) {
                            if (d[j] < e0) {
                                e0 = d[j];
                                col = j;
                            }
                        }

                        System.out.println("Iteri = " + iter + "  Ener = " + e0);
                        //if this is the end of iteration if Energy (e1-e0)<tolerance (1.e-7)....

                        if ((force_maxiter && iter >= control_max) || (iter >= dim - 1 || iter == 999 || Math.abs(b[iter]) < tolerance) || (!force_maxiter && Math.abs(eini - e0) <= tolerance)) {
                            // converged
                            System.out.println("E0 = " + e0);
                            maxiter = iter;
                            if (!calc_gs) {
                                return; // We return with ground states energy
                            }
                            z = Array.Identity(z.length); //identity;
                            Array.copy(d, 1, iter, a, 1, iter);
                            Array.copy(e, 2, iter + 1, b, 1, iter);
                            // call tqli with eigenvectors
                            tqli(d, e, iter, z, true);
                            xc = new ArrayRealVector(Array.copy(z[col], 0, iter));

                            Array.copy(d, 0, maxiter - 1, d, 1, maxiter);
                            printEnergyLevel(maxiter, d);

                            gs.set(0.0);
                            ener = e0;
                            for (int i = 0; i < maxiter; i++) {
                                x2 = save.get(i);
                                gs = gs.add(x2.mapMultiply(xc.getEntry(i)));
                            }

                            return;
                        }
                    } // diagonalization of tridiagonal matrix
                }
            } // Lanczos iteration
        } // main iteration

    }

    /* (C) Copr. 1986-92 Numerical Recipes Software . */
    private void tqli(double[] d, double[] e, int n, double[][] eigVec, boolean calc_eigVec) {
        int m, l, iter, i, k;
        double s, r, p, g, f, dd, c, b;

        for (i = 2; i <= n; i++) {
            e[i - 1] = e[i];
        }
        e[n] = 0.0f;
        for (l = 1; l <= n; l++) {
            iter = 0;
            do {
                for (m = l; m <= n - 1; m++) {
                    dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
                    if (Math.abs(e[m]) + dd == dd)
                        break;
                }
                if (m != l) {
                    if (iter++ == 30) {
                        System.out.println("Too many iterations in TQLI");
                    }
                    g = (d[l + 1] - d[l]) / (2.0f * e[l]);
                    r = Math.sqrt((g * g) + 1.0f);
                    g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                    s = c = 1.0f;
                    p = 0.0f;
                    for (i = m - 1; i >= l; i--) {
                        f = s * e[i];
                        b = c * e[i];
                        if (Math.abs(f) >= Math.abs(g)) {
                            c = g / f;
                            r = Math.sqrt((c * c) + 1.0f);
                            e[i + 1] = f * r;
                            c *= (s = 1.0f / r);
                        } else {
                            s = f / g;
                            r = Math.sqrt((s * s) + 1.0f);
                            e[i + 1] = g * r;
                            s *= (c = 1.0f / r);
                        }
                        g = d[i + 1] - p;
                        r = (d[i] - g) * s + 2.0f * c * b;
                        p = s * r;
                        d[i + 1] = g + p;
                        g = c * r - b;
                        /* Next loop can be omitted if eigenvectors not wanted */
                        if (calc_eigVec) {
                            for (k = 1; k <= n; k++) {
                                f = eigVec[k][i + 1];
                                eigVec[k][i + 1] = s * eigVec[k][i] + c * f;
                                eigVec[k][i] = c * eigVec[k][i] - s * f;
                            }
                        }
                    }
                    d[l] = d[l] - p;
                    e[l] = g;
                    e[m] = 0.0f;
                }
            } while (m != l);
        }
    }

    private void printEnergyLevel(int maxIter, double[] alpha) {
        double[] copy = alpha;
        double temp;
        if (maxIter == 1) {
            System.out.println(alpha[0]);
        }
        for (int i = 0; i < maxIter; i++) {
            for (int j = i + 1; j < maxIter; j++) {
                if (copy[i] > copy[j]) {
                    temp = copy[i];
                    copy[i] = copy[j];
                    copy[j] = temp;
                }
            }
            System.out.println("Energy LEVEL " + i + " : " + copy[i]);
        }
        System.out.println("===================================================");
    }
}
