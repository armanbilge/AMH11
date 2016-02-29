/**
 * Utils.java
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
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;

import java.util.Random;

/**
 * @author Arman Bilge
 */
public final class Utils {

    private Utils() {}

    public interface MatrixFunction {
        Matrix apply(Matrix X, boolean transpose);
        int getDimensions();
        boolean isReal();
    }

    public interface AMH11Random {
        boolean nextBoolean();
    }

    private static AMH11Random random = new AMH11Random() {
        private final Random random = new Random();
        public boolean nextBoolean() {
            return random.nextBoolean();
        }
    };

    public static void setRandom(AMH11Random random) {
        Utils.random = random;
    }

    public static Vector fill(Vector v, double d) {
        for (int i = 0; i < v.size(); ++i) v.set(i, d);
        return v;
    }

    public static double trace(Matrix M) {
        double trace = 0.0;
        for (int i = 0; i < M.numRows(); ++i)
            trace += M.get(i, i);
        return trace;
    }

    public static Matrix ascendingDiagonal(int n) {
        Matrix M = new DenseMatrix(n, n);
        for (int i = 0; i < n; ++i) M.set(i, i, i+1);
        return M;
    }

    public static Matrix ceil(Matrix M) {
        for (MatrixEntry me : M) me.set(Math.ceil(me.get()));
        return M;
    }

    public static Matrix abs(Matrix M) {
        for (MatrixEntry me : M) me.set(Math.abs(me.get()));
        return M;
    }

    public static Matrix Zero2Inf(Matrix M) {
        for (int i = 0; i < M.numRows(); ++i) {
            for (int j = 0; j < M.numColumns(); ++j)
                if (M.get(i, j) == 0) M.set(i, j, Double.POSITIVE_INFINITY);
        }
        return M;
    }

    public static boolean isPositive(Matrix M) {
        for (MatrixEntry me : M) {
            if (me.get() < 0) return false;
        }
        return true;
    }

    public static double entrySum(Vector v) {
        double sum = 0.0;
        for (VectorEntry ve : v) sum += ve.get();
        return sum;
    }

    public static Vector rowAsVector(Matrix M, int row) {
        Vector v = new DenseVector(M.numColumns());
        for (int i = 0; i < M.numColumns(); ++i) v.set(i, M.get(row, i));
        return v;
    }

    public static Vector columnAsVector(Matrix M, int column) {
        Vector v = new DenseVector(M.numRows());
        for (int i = 0; i < M.numColumns(); ++i) v.set(i, M.get(i, column));
        return v;
    }

    public static Vector absColumnSums(Matrix M) {
        Vector sums = new DenseVector(M.numColumns());
        for (int i = 0; i < M.numColumns(); ++i) {
            double sum = Utils.columnAsVector(M, i).norm(Vector.Norm.One);
            sums.set(i, sum);
        }
        return sums;
    }

    public static Matrix randomSigns(int r, int c) {
        Matrix signs = new DenseMatrix(r, c);
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                double sign = random.nextBoolean() ? 1.0 : -1.0;
                signs.set(i, j, sign);
            }
        }
        return signs;
    }

    public static Vector columnsMax(Matrix M) {
        Vector maxes = new DenseVector(M.numColumns());
        for (int i = 0; i < maxes.size(); ++i)
            maxes.set(i, Utils.columnAsVector(M, i).norm(Vector.Norm.Infinity));
        return maxes;
    }

    public static int memberCount(int[] a, Vector b) {
        int count = 0;
        for (int d : a)
            if (contains(b, d)) ++count;
        return count;
    }

    public static boolean contains(Vector v, double d) {
        for (int i = 0; i < v.size(); ++i) {
            if (v.get(i) == d) return true;
        }
        return false;
    }

    public static Matrix createSubMatrix(Matrix M, int row, int column,
            int rows, int columns) {
        Matrix S = new DenseMatrix(rows, columns);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j)
                S.set(i, j, M.get(row+i, column+j));
        }
        return S;
    }
}
