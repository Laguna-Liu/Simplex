package com.company;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.Algebra;

public class SimplexMethod {

    public void print(Object T){
        System.out.println(T);
    }

    public DoubleMatrix2D m_multi(DoubleMatrix2D left, DoubleMatrix2D right){
        return Algebra.DEFAULT.mult(left, right);
    }

    public DoubleMatrix2D solve(DoubleMatrix2D A, DoubleMatrix2D b, DoubleMatrix2D c){
        // initialize

        DoubleFactory2D factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D x = new DenseDoubleMatrix2D(A.columns(), 1);

        DoubleMatrix2D B = factory2D.identity(A.rows());
        DoubleMatrix2D c_b = new DenseDoubleMatrix2D(A.rows(), 1);

        int iterations = 0;

        while(iterations++ < 5000){

            print("Iteration = " + iterations);

            DoubleMatrix2D c_b_t = Algebra.DEFAULT.transpose(c_b);
            DoubleMatrix2D B_i = Algebra.DEFAULT.inverse(B);
            DoubleMatrix2D w_t = Algebra.DEFAULT.mult(c_b_t, B_i);
            DoubleMatrix2D x_b = m_multi(B_i, b);

            double[][] z_values = new double[x.rows()][1];

            for(int j = 0; j < x.rows(); j ++) {
                z_values[j][0] = m_multi(w_t, A.viewPart(0, j, A.columns(),1)).toArray()[0][0];
            }

            DoubleMatrix2D z = new DenseDoubleMatrix2D(z_values.length, 1);
            DoubleMatrix2D z0 = m_multi(w_t, b);
            z.assign(z_values);

            double max_z_c = z.toArray()[0][0] - c.toArray()[0][0];
            int k = 0;
            for(int i = 0; i < x.rows(); i ++){
                double z_c = z.toArray()[i][0] - c.toArray()[i][0];
                if(z_c > max_z_c) {
                    max_z_c = z_c;
                    k = i;
                }
            }

            double[] b_s = new double[x_b.rows()];
            for(int i = 0; i < x_b.rows(); i ++){
                b_s[i] = x_b.toArray()[i][0];
            }
            print(x_b);
            DoubleMatrix2D yk = m_multi(B_i, A.viewPart(0, k, A.columns(),1));
            double[] ys_k = new double[yk.rows()];
            for(int i = 0; i < ys_k.length; i ++)
                ys_k[i] = yk.toArray()[i][0];

            double x_k = 999999;
            int r = -1;

            for(int i = 0; i < x_b.rows(); i ++){
                if(ys_k[i] > 0 && b_s[i] / ys_k[i] < x_k){
                    r = i;
                    x_k = b_s[i] / ys_k[i];
                }
            }

            for(int i = 0; i < A.rows(); i ++){
                double ai = A.get(i, k);
                double bi = B.get(i, r);
                A.set(i, k, bi);
                B.set(i, r, ai);
            }

            double ci = c.get(k, 0);
            double c_b_i = c_b.get(r, 0);

            c.set(k, 0, c_b_i);
            c_b.set(r, 0, ci);

            if(max_z_c <= 0){
                print("find maximum");
                print(z0);
                return z0;
            }

        }


        return null;
    }

    public static void insert(){

        SimplexMethod simplexMethod = new SimplexMethod();

        DoubleMatrix2D A = new DenseDoubleMatrix2D(2, 2);
        DoubleMatrix2D b = new DenseDoubleMatrix2D(2, 1);
        DoubleMatrix2D c = new DenseDoubleMatrix2D(2, 1);

        A.assign(new double[][]{{-1, 1}, {1, 1}});
        b.assign(new double[][]{{1}, {4}});
        c.assign(new double[][]{{-1}, {-2}});

        simplexMethod.solve(A, b, c);


    }

    public static void main(String[] args){

//        SimplexMethod simplexMethod = new SimplexMethod();
//
//        DoubleMatrix2D A = new DenseDoubleMatrix2D(3, 3);
//        DoubleMatrix2D b = new DenseDoubleMatrix2D(3, 1);
//        DoubleMatrix2D c = new DenseDoubleMatrix2D(3, 1);
//
//        A.assign(new double[][]{{3, -1, 2}, {-2, 4, 0}, {-4, 3, 8}});
//        b.assign(new double[][]{{7}, {12}, {10}});
//        c.assign(new double[][]{{1}, {-3}, {2}});
//
//        simplexMethod.solve(A, b, c);

        insert();

    }


}
