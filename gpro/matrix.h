#ifndef MATRIX_H
#define MATRIX_H
#include <iomanip>
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
// Various helper routines for debugging
// Output a matrix
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& matrix);

// Matrix multiplication
std::vector<std::vector<double>> MatrixTimesMatrix(const std::vector<std::vector<double>>& A,
                                                   const std::vector<std::vector<double>>& B);
// Matrix times vector
std::vector<double> MatrixTimesVector(const std::vector<std::vector<double>>& A, const std::vector<double>& x);
// Cholesky matrix decomposition to lower triangular matrix and its conjugate transpose
// 
// Restricted to positive-definite matrices
class CholeskyDecomposition {
public:
    // Constructor
    // Matrix is decomposed in-place
    CholeskyDecomposition(std::vector<std::vector<double>>& sourceMatrix);

    // Destructor
    ~CholeskyDecomposition();

    // Decomposition into triangular matrices
    bool Decompose();

    // Solve for x in form Ax = b.  A is the original input matrix.
    std::vector<double> Solve(const std::vector<double>& b);

protected:
    // Input matrix
    std::vector<std::vector<double>>& decomposedMatrix;

private:

    CholeskyDecomposition(const CholeskyDecomposition&);

    void operator=(const CholeskyDecomposition&);
};

class LUDecomposition {
public:
    // Constructor
    // Matrix is decomposed in-place
    LUDecomposition(std::vector<std::vector<double>>& sourceMatrix);

    // Destructor
    ~LUDecomposition();

    // Decomposition into triangular matrices
    bool Decompose();

    // Solve for x in form Ax = b.  A is the original input matrix.
    std::vector<double> Solve(const std::vector<double>& b);

protected:
    // Output matrix after decomposition
    std::vector<std::vector<double>>& decomposedMatrix;

    // Permutation of rows during pivoting
    std::vector<int> rowPermutation;

private:
    LUDecomposition(const LUDecomposition&);
    void operator=(const LUDecomposition&);
};
#endif
