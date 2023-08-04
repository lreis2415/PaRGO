#include "matrix.h"
// Various helper routines for debugging
// Output a matrix
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& matrix) {
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix.size(); j++) {
            if (j != 0) {
                std::cout << " ";
            }
            std::cout << std::setw(10) << std::setprecision(7) << matrix[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return os;
}

// Matrix multiplication
std::vector<std::vector<double>> MatrixTimesMatrix(const std::vector<std::vector<double>>& A,
                                                   const std::vector<std::vector<double>>& B) {
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < B[0].size(); j++) {
            for (size_t k = 0; k < A[0].size(); k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Matrix times vector
std::vector<double> MatrixTimesVector(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    std::vector<double> b(x.size(), 0.0);
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < x.size(); j++) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}

CholeskyDecomposition::CholeskyDecomposition(std::vector<std::vector<double>>& sourceMatrix) : decomposedMatrix(sourceMatrix) {
    assert(sourceMatrix.size() > 0 && sourceMatrix.size() == sourceMatrix[0].size());
}

// Destructor
CholeskyDecomposition::~CholeskyDecomposition() {
}

// Decomposition into triangular matrices
bool CholeskyDecomposition::Decompose() {
    // Enumerate matrix columnwise
    for (size_t j = 0; j < decomposedMatrix.size(); j++) {
        for (size_t i = j; i < decomposedMatrix.size(); i++) {
            if (i == j) {
                double sum = 0.0;
                for (size_t k = 0; k < i; k++) {
                    sum += std::pow(decomposedMatrix[i][k], 2.0);
                }
                if (decomposedMatrix[i][j] - sum <= 0.0) {
                    // Not positive definite matrix
                    return false;
                }
                decomposedMatrix[i][j] = std::sqrt(decomposedMatrix[i][j] - sum);
            }
            else {
                double sum = 0.0;
                for (size_t k = 0; k < j; k++) {
                    sum += (decomposedMatrix[i][k] * decomposedMatrix[j][k]);
                }
                decomposedMatrix[i][j] = (1 / decomposedMatrix[j][j]) * (decomposedMatrix[i][j] - sum);
                decomposedMatrix[j][i] = decomposedMatrix[i][j];
            }
        }
    }
    return true;
}

// Solve for x in form Ax = b.  A is the original input matrix.
std::vector<double> CholeskyDecomposition::Solve(const std::vector<double>& b) {
    std::vector<double> y(b.size());
    // First solve lower triangular * y = b with forward substitution
    for (size_t i = 0; i < b.size(); i++) {
        double sum = 0.0;
        for (size_t j = 0; j < i; j++) {
            sum += (decomposedMatrix[i][j] * y[j]);
        }
        y[i] = (b[i] - sum) / decomposedMatrix[i][i];
    }
    // Now solve upper triangular (transpose of lower triangular) * x = y with back substitution.
    // Note that x can be solved in place using the existing y vector.  No need to allocate 
    // another vector.
    for (int i = static_cast<int>(b.size()) - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = static_cast<int>(b.size()) - 1; j > i; j--) {
            sum += (decomposedMatrix[i][j] * y[j]);
        }
        y[i] = (y[i] - sum) / decomposedMatrix[i][i];
    }
    return y;
}

// Constructor
LUDecomposition::LUDecomposition(std::vector<std::vector<double>>& sourceMatrix) : decomposedMatrix(sourceMatrix) {
    assert(sourceMatrix.size() > 0 && sourceMatrix.size() == sourceMatrix[0].size());
}

// Destructor
LUDecomposition::~LUDecomposition() {
}

// Decomposition into triangular matrices
bool LUDecomposition::Decompose() {
    // Initialize the permutation vector
    size_t n = decomposedMatrix.size();
    rowPermutation.reserve(n);
    for (size_t i = 0; i < n; i++) {
        rowPermutation.push_back(static_cast<int>(i));
    }
    double det = 1.0;
    // LU factorization.
    for (size_t p = 1; p <= n - 1; p++) {
        // Find pivot element.
        for (size_t i = p + 1; i <= n; i++) {
            if (std::fabs(decomposedMatrix[rowPermutation[i - 1]][p - 1]) > std::fabs(decomposedMatrix[rowPermutation[p - 1]][p - 1])) {
                // Switch the index for the p-1 pivot row if necessary.
                std::swap(rowPermutation[p - 1], rowPermutation[i - 1]);
                det = -det;
            }
        }
        if (decomposedMatrix[rowPermutation[p - 1]][p - 1] == 0.0) {
            // The matrix is singular, at least to precision of algorithm
            return false;
        }
        // Multiply the diagonal elements.
        det = det * decomposedMatrix[rowPermutation[p - 1]][p - 1];
        // Form multiplier.
        for (size_t i = p + 1; i <= n; i++) {
            decomposedMatrix[rowPermutation[i - 1]][p - 1] /= decomposedMatrix[rowPermutation[p - 1]][p - 1];
            // Eliminate [p-1].
            for (size_t j = p + 1; j <= n; j++) {
                decomposedMatrix[rowPermutation[i - 1]][j - 1] -= decomposedMatrix[rowPermutation[i - 1]][p - 1] * decomposedMatrix[rowPermutation[p - 1]][j - 1];
            }
        }
    }
    det = det * decomposedMatrix[rowPermutation[n - 1]][n - 1];
    return (det != 0.0);
}

// Solve for x in form Ax = b.  A is the original input matrix.
// Note: b is modified in-place for row permutations
std::vector<double> LUDecomposition::Solve(const std::vector<double>& b) {
    // Our decomposed matrix is comprised of both the lower and upper diagonal matrices.
    // The rows of this matrix have been permutated during the decomposition process.  The
    // rowPermutation indicates the proper row order.
    // The lower diagonal matrix only include elements below the diagonal with diagonal 
    // elements set to 1.
    // The upper diagonal matrix is fully specified.
    // First solve Ly = Pb for y using forward substitution. P is a permutated identity matrix.
    std::vector<double> y(b.size());
    for (size_t i = 0; i < y.size(); i++) {
        size_t currentRow = rowPermutation[i];
        double sum = 0.0;
        for (size_t j = 0; j < i; j++) {
            sum += (decomposedMatrix[currentRow][j] * y[j]);
        }
        y[i] = (b[currentRow] - sum);
    }
    // Now solve Ux = y for x using back substitution.  Note that 
    // x can be solved in place using the existing y vector.  No need
    // to allocate another vector.
    for (int i = static_cast<int>(b.size()) - 1; i >= 0; i--) {
        size_t currentRow = rowPermutation[i];
        double sum = 0.0;
        for (int j = static_cast<int>(b.size()) - 1; j > i; j--) {
            sum += (decomposedMatrix[currentRow][j] * y[j]);
        }
        y[i] = (y[i] - sum) / decomposedMatrix[currentRow][i];
    }
    return y;
}
