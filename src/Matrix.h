#pragma once
#include <cstddef>
#include <iosfwd>
#include <vector>

// Плотная матрица вещественных чисел и базовые операции над ней.
class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, double value = 0.0);

    static Matrix identity(std::size_t n);

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }

    double& operator()(std::size_t r, std::size_t c);
    double operator()(std::size_t r, std::size_t c) const;

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;

    Matrix transposed() const;
    double determinant() const;   // метод Гаусса с выбором главного элемента
    Matrix inverse() const;       // метод Гаусса–Жордана

private:
    std::size_t rows_ = 0, cols_ = 0;
    std::vector<double> data_;
};

std::ostream& operator<<(std::ostream& os, const Matrix& m);
