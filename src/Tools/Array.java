package Tools;

import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;

public class Array {   
    
    public static double[] getRandom(int d) {
        UnitSphereRandomVectorGenerator uvg = new UnitSphereRandomVectorGenerator(d);
        return uvg.nextVector();
    }
    
    public static void setValues(double[] arr, double val) {
        for(int i=0; i<arr.length; i++) arr[i] = val;
    }
    
    
    public static double[] setValue(double[] a, double[] b, int begin, int end) {
        if((end-begin) > b.length) {
            System.err.println("range index array a to short");
            return null;
        }
        for(int i=0; i<b.length; i++)
            a[begin+i] = b[i];
        return a;
    }
    
    /**
     * Copy value from array B in range (Bo,Bt) to array A in range (Ao,At) ;
     * @param A the main array;
     * @param B the source of values
     * @param Ao first index array A
     * @param At the last index array A
     */
    public static void copy(double[] A, int Ao, int At, double[] B, int Bo, int Bt) {
        int lengA = At-Ao;
        int lengB = Bt-Bo;
        if(lengA != lengB) {System.err.println("range index array doesn't match");} 
        else
            for(int i=0; i<=lengB; i++)
                A[Ao+i] = B[Bo+i];
    }
    
    /**
     * Copy value from array B in range (Bo,Bt) to all elements in array A;
     * @param a the source of values
     * @param start first index array B
     * @param end the last index array B
     */
    public static double[] copy(double[] a, int start, int end) {
        int lenght = (end-start)+1;
        double[] x = new double[lenght];
        for(int i=0; i<lenght; i++)
            x[i] = a[start+i];
        return x;
    }
    
    public static double[][] Identity(int dim) {
        double[][] I = new double[dim][dim];
        for(int i=0; i<dim; i++) I[i][i] = 1.0;
        return I;
    }
    
    public static double[] plus(double[] a, double[] b) {
        if(!isValid(a, b)) return null;
        double[] result = new double[a.length];
        for(int i=0; i<a.length; i++) result[i] = a[i]+b[i];
        return result;
    }
    
    public static double[] minus(double[] a, double[] b) {
        if(!isValid(a, b)) return null;
        double[] result = new double[a.length];
        for(int i=0; i<a.length; i++) result[i] = a[i]-b[i];
        return result;
    }
    
    /**************************** PRODUCT *************************************/
    public static double getNorm(double[] a) {return product(a,a);}
    
    /**
     * [scalar product] multiply each element of arr with scalar x
     * @param arr an array double[]
     * @param x double scalar
     * @return new array double[]
     */
    public static double[] product(double[] arr, double x) {
        double[] result = new double[arr.length];
        for(int i=0; i<arr.length; i++) result[i] = arr[i] * x;
        return result;
    }
    
    /**
     * [dot product] array a and b
     * @param a array
     * @param b array
     * @return new double
     */
    public static double product(double[] a, double[] b) {
        double result = 0;
        for(int i=0; i<a.length; i++) result += a[i]*b[i];
        return result;
    }
    
    public static double[] product(double[][] a, double[] b) {
        if(a[0].length != b.length) {
            System.err.println("Matrix/Array lenght doesnt match");
            return null;
        }
        
        double[] result = new double[a.length];
        for(int i=0; i<a.length; i++) 
            for(int j=0; j<a[0].length; j++)
                result[i] += a[i][j]*b[j];
        return result;
    }
    
    /***************************  DIVIDE **************************************/
    /**
     * [scalar product] divide each element of arr with scalar x
     * @param arr an array double[]
     * @param x double scalar
     * @return new array double[]
     */
    public static double[] divide(double[] arr, double x) {
        double[] result = new double[arr.length];
        for(int i=0; i<arr.length; i++) result[i] = arr[i] / x;
        return result;
    }
       
    public static double[] shortDesc(double[] v) {
        double temp;
        for(int i=0; i<v.length-1; i++) {
            temp = v[i];
            for(int j=i+1; j<v.length; j++) {
                if(temp > v[j]) {
                    v[i] = v[j];
                    v[j] = temp;
                    temp = v[i];
                }
            }
        }
        return v;
    }
    
    private static boolean isValid(double[] a, double[] b) {
        if(a.length != b.length) {
            System.err.println("lenght array douesn't match");
            return false;
        }
        else
            return true;
    }
    
    public static void print(double[] v) {
        for (int i = 0; i < v.length; i++) 
            System.out.println(v[i]);
    }
}
