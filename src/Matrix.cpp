#include "Matrix.h"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace {
constexpr double kEps = 1e-12;
}

Matrix::Matrix(std::size_t rows, std::size_t cols, double value)
    : rows_(rows), cols_(cols), data_(rows * cols, value) {}

Matrix Matrix::identity(std::size_t n) {
    Matrix m(n, n);
    for (std::size_t i = 0; i < n; ++i) m(i, i) = 1.0;
    return m;
}

double& Matrix::operator()(std::size_t r, std::size_t c) {
    if (r >= rows_ || c >= cols_) throw std::out_of_range("Matrix index out of range");
    return data_[r * cols_ + c];
}

double Matrix::operator()(std::size_t r, std::size_t c) const {
    if (r >= rows_ || c >= cols_) throw std::out_of_range("Matrix index out of range");
    return data_[r * cols_ + c];
}

Matrix Matrix::operator+(const Matrix& o) const {
    if (rows_ != o.rows_ || cols_ != o.cols_) throw std::invalid_argument("Sizes differ");
    Matrix res(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] + o.data_[i];
    return res;
}

Matrix Matrix::operator-(const Matrix& o) const {
    if (rows_ != o.rows_ || cols_ != o.cols_) throw std::invalid_argument("Sizes differ");
    Matrix res(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] - o.data_[i];
    return res;
}

Matrix Matrix::operator*(const Matrix& o) const {
    if (cols_ != o.rows_) throw std::invalid_argument("Incompatible sizes for multiplication");
    Matrix res(rows_, o.cols_);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t k = 0; k < cols_; ++k) {
            const double a = (*this)(i, k);
            for (std::size_t j = 0; j < o.cols_; ++j) res(i, j) += a * o(k, j);
        }
    return res;
}

Matrix Matrix::operator*(double s) const {
    Matrix res(*this);
    for (double& x : res.data_) x *= s;
    return res;
}

Matrix Matrix::transposed() const {
    Matrix res(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t j = 0; j < cols_; ++j) res(j, i) = (*this)(i, j);
    return res;
}

double Matrix::determinant() const {
    if (rows_ != cols_) throw std::invalid_argument("Determinant requires a square matrix");
    Matrix a(*this);
    double det = 1.0;
    for (std::size_t col = 0; col < cols_; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < rows_; ++r)
            if (std::fabs(a(r, col)) > std::fabs(a(pivot, col))) pivot = r;
        if (std::fabs(a(pivot, col)) < kEps) return 0.0;
        if (pivot != col) {
            for (std::size_t j = 0; j < cols_; ++j) std::swap(a(pivot, j), a(col, j));
            det = -det;
        }
        det *= a(col, col);
        for (std::size_t r = col + 1; r < rows_; ++r) {
            const double f = a(r, col) / a(col, col);
            for (std::size_t j = col; j < cols_; ++j) a(r, j) -= f * a(col, j);
        }
    }
    return det;
}

Matrix Matrix::inverse() const {
    if (rows_ != cols_) throw std::invalid_argument("Inverse requires a square matrix");
    const std::size_t n = rows_;
    Matrix a(*this), inv = identity(n);
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::fabs(a(r, col)) > std::fabs(a(pivot, col))) pivot = r;
        if (std::fabs(a(pivot, col)) < kEps) throw std::domain_error("Matrix is singular");
        for (std::size_t j = 0; j < n; ++j) {
            std::swap(a(pivot, j), a(col, j));
            std::swap(inv(pivot, j), inv(col, j));
        }
        const double p = a(col, col);
        for (std::size_t j = 0; j < n; ++j) { a(col, j) /= p; inv(col, j) /= p; }
        for (std::size_t r = 0; r < n; ++r) {
            if (r == col) continue;
            const double f = a(r, col);
            for (std::size_t j = 0; j < n; ++j) {
                a(r, j) -= f * a(col, j);
                inv(r, j) -= f * inv(col, j);
            }
        }
    }
    return inv;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) os << std::setw(10) << std::setprecision(4) << m(i, j);
        os << '\n';
    }
    return os;
}
