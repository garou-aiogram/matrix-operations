// Временная консольная проверка ядра. Оконный интерфейс (Windows) будет добавлен позже.
#include <iostream>
#include "Matrix.h"

int main() {
    Matrix a(2, 2);
    a(0, 0) = 4; a(0, 1) = 7;
    a(1, 0) = 2; a(1, 1) = 6;

    std::cout << "A:\n" << a
              << "det(A) = " << a.determinant() << "\n"
              << "A^T:\n" << a.transposed()
              << "A^-1:\n" << a.inverse()
              << "A * A^-1:\n" << a * a.inverse();
}
