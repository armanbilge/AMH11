/**
 * HT00.java
 *
 * AMH11: Java implementation of the matrix exponential method
 *     described by Al-Mohy and Higham (2011)
 *
 * Copyright (c) 2016 Arman Bilge
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.compevol.amh11;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import org.compevol.amh11.Utils.MatrixFunction;

import java.util.Arrays;
import java.util.Comparator;

/**
 * @author Arman Bilge
 */
public final class HT00 {

    private HT00() {}

    public static double[] normest1(MatrixFunction f, int t) {

        int n = f.getDimensions();
        boolean isReal = f.isReal();

        t = Math.abs(t);

//        int rpt_S = 0;
//        int rpt_e = 0;

        if (t == n || n <= 4) {
            Matrix X = Matrices.identity(n);
            Matrix Y = f.apply(X, false);

            Vector sums = Utils.absColumnSums(Y);
            double est = sums.norm(Vector.Norm.Infinity);
            return new double[]{est, 1};
        }

        Matrix X = Utils.randomSigns(n, t);
        for (int i = 0; i < n; ++i) X.set(i, 0, 1.0);
        unduplicate(X, null);
        X.scale(1.0 / n);

        int itMax = 5;
        int it = 0;
        int nmv = 0;

        Vector ind = new DenseVector(t);
        Vector vals;
        Vector Zvals = new DenseVector(n);
        Matrix S = new DenseMatrix(n, t);
        Matrix S_old;
        Vector indHist = null;

        double est;
        double estOld = 0.0;

        while (true) {

            ++it;
            Matrix Y = f.apply(X, false);
            ++nmv;

            vals = Utils.absColumnSums(Y);
            final double[] valsArray = Matrices.getArray(vals);
            Integer[] m = new Integer[vals.size()];
            for (int i = 0; i < m.length; ++i) m[i] = i;
            Arrays.sort(m, new Comparator<Integer>() {
                public int compare(Integer i, Integer j) {
                    return Double.compare(valsArray[i], valsArray[j]);
                }
            });
            int[] m2 = new int[t];
            for (int i = 0; i < t; ++i) {
                m2[i] = m[t-1-i];
            }
            vals = Matrices.getSubVector(vals, m2);
            Vector vals_ind = Matrices.getSubVector(ind, m2);
            est = vals.get(0);

            int est_j = Integer.MIN_VALUE;
            if (est > estOld || (it == 2)) {
                est_j = (int) vals_ind.get(0);
            }

            if (it >= 2 && est <= estOld) {
                est = estOld;
                break;
            }
            estOld = est;

            if (it > itMax) {
//                it = itMax;
                break;
            }

            S_old = S;
            S = mySign(Y);

            if (isReal) {
                Matrix SS = S_old.transAmult(S,
                        new DenseMatrix(S_old.numColumns(), S.numColumns()));
                double np = Utils.entrySum(Utils.columnsMax(SS));
                if (np == t) break;
                /* double r = */ unduplicate(S, S_old);
//                rpt_S += r;
            }

            Matrix Z = f.apply(S, true);
            ++nmv;
            for (int i = 0; i < n; ++i)
                Zvals.set(i,
                        Utils.rowAsVector(Z, i).norm(Vector.Norm.Infinity));

            if (it >= 2 && Zvals.norm(Vector.Norm.Infinity) == Zvals.get(est_j))
                break;

            m = new Integer[Zvals.size()];
            final double[] zvalsArray = Matrices.getArray(Zvals);
            for (int i = 0; i < m.length; ++i) m[i] = i;
            Arrays.sort(m, new Comparator<Integer>() {
                public int compare(Integer i, Integer j) {
                    return Double.compare(zvalsArray[i], zvalsArray[j]);
                }
            });
            m2 = new int[n];
            for (int i = 0; i < n; ++i) m2[i] = m[n-1-i];

            int imax = t;
            int[] m2Tot = Arrays.copyOfRange(m2, 0, t);
            if (it == 1) {
                ind = new DenseVector(t);
                for (int i = 0; i < t; ++i) ind.set(i, m2Tot[i]);
                indHist = ind;
            } else {
                int rep = Utils.memberCount(m2Tot, indHist);
//                rpt_e += rep;
                if (rep == t) break;
                int j = 0;
                for (int i = 0; i < t; ++i) {
                    if (j >= n) {
                        imax = i-1;
                        break;
                    }
                    while (Utils.contains(indHist, m2[j])) {
                        ++j;
                        if (j >= n) {
                            imax = i - 1;
                            break;
                        }
                    }
                    if (j >= n) break;
                    ind.set(i, m2[j]);
                    ++j;
                }
                Vector temp = new DenseVector(t + imax);
                for (int i = 0; i < t; ++i)
                    temp.set(i, indHist.get(i));
                for (int i = 0; i < imax; ++i)
                    temp.set(t+i, ind.get(i));
                indHist = temp;
            }
            X = new DenseMatrix(n, t);
            for (int j = 0; j < imax; ++j)
                X.set((int) ind.get(j), j, 1);
        }

        int iter = nmv;
        return new double[]{est, iter};

    }

    private static double unduplicate(Matrix S,
            Matrix S_old) {

        int n = S.numRows();
        int t = S.numColumns();
        double r = 0.0;
        if (t == 1) return r;

        Matrix W;
        int jStart;
        int lastCol;
        if (S_old == null) {
            W = new DenseMatrix(n, t);
            for (int i = 0; i < n; ++i)
                W.set(i, 0, S.get(i, 0));
            jStart = 1;
            lastCol = 1;
        } else {
            W = new DenseMatrix(n, 2*t-1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < t; ++j) {
                    W.set(i, j, S_old.get(i, j));
                }
            }
            jStart = 0;
            lastCol = t;
        }

        for (int j = jStart; j < t; ++j) {
            int rpt = 0;
            while (Utils.columnAsVector(Utils.abs(
                    Utils.createSubMatrix(S, 0, j, n, 1)
                    .transAmult(Utils.createSubMatrix(W, 0, 0, n, lastCol),
                            new DenseMatrix(1, lastCol))), 0)
                            .norm(Vector.Norm.Infinity) == n) {
                ++rpt;
                Matrix signs = Utils.randomSigns(n, 1);
                for (int i = 0; i < n; ++ i)
                    S.set(i, j, signs.get(i, 0));
                if (rpt > n / (double) t) break;
            }

            r += Math.signum(rpt);
            if (j < t) {
                ++lastCol;
                for (int i = 0; i < n; ++i) {
                    W.set(i, lastCol-1, S.get(i, j));
                }
            }
        }

        return r;
    }

    private static Matrix mySign(Matrix M) {
        Matrix signs = new DenseMatrix(M.numRows(), M.numColumns());
        for (int i = 0; i < M.numRows(); ++i) {
            for (int j = 0; j < M.numColumns(); ++j) {
                double sign = Math.signum(M.get(i, j));
                if (sign == 0) sign = 1.0;
                signs.set(i, j, sign);
            }
        }
        return signs;
    }

}
