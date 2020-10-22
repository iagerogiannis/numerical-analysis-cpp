#include "Polynomial.h"

#include <iostream>
#include <cmath>
#include "cpp_extended.h"
#include "math_extended.h"

Polynomial::Polynomial() = default;

Polynomial::Polynomial(int n, double* coefficients_)
    : n(n) {
    coefficients = copyDynamicArray(coefficients_, n + 1);
}

Polynomial::Polynomial(const Polynomial & polynomial) {

    n = polynomial.n;
    coefficients = copyDynamicArray(polynomial.coefficients, n + 1);
}

Polynomial::Polynomial(Polynomial && polynomial) noexcept {

    n = polynomial.n;
    coefficients = polynomial.coefficients;

    polynomial.coefficients = nullptr;
}

Polynomial& Polynomial::operator=(Polynomial && polynomial)  noexcept {

    if (this != &polynomial) {
        delete[] coefficients;

        n = polynomial.n;
        coefficients = polynomial.coefficients;

        polynomial.coefficients = nullptr;
    }
    return *this;
}

Polynomial::~Polynomial() {
    delete[] coefficients;
    coefficients = nullptr;
}

void Polynomial::setCoefficients(int n_, double* coefficients_) {
    n = n_;
    coefficients = copyDynamicArray(coefficients_, n + 1);
}

double Polynomial::value(double x) {
    double result = 0;
    for (int i = 0; i < n + 1; ++i) {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

double Polynomial::definite_integral(double a, double b) {
    double integral_value = 0.;
    for (int i = 0; i < n + 1; ++i) {
        integral_value += (coefficients[i] / (double)(i + 1)) * (pow(b, i + 1) - pow(a, i + 1));
    }
    return integral_value;
}

Polynomial Polynomial::derivative() {
    auto* derivative_coefficients = new double[n];
    for (int i = 1; i < n + 1; ++i) {
        derivative_coefficients[i - 1] = i * coefficients[i];
    }
    return { n - 1, derivative_coefficients };
}

Polynomial Polynomial::indefinite_integral() {
    auto* integral_coefficients = new double[n + 2];
    integral_coefficients[0] = 0.;
    for (int i = 0; i < n + 1; ++i) {
        integral_coefficients[i + 1] = coefficients[i] / (double)(i + 1);
    }
    return { n + 1, integral_coefficients };
}

Polynomial Polynomial::operator+(Polynomial polynomial) {

    int big_n = std::max(n, polynomial.n);
    int small_n = std::min(n, polynomial.n);

    double* result_coefficients;
    result_coefficients = new double[big_n];

    for (int i = 0; i < small_n + 1; ++i) {
        result_coefficients[i] = coefficients[i] + polynomial.coefficients[i];
    }

    if (n > small_n) {
        for (int i = small_n + 1; i < big_n + 1; ++i) {
            result_coefficients[i] = coefficients[i];
        }
    }
    else {
        for (int i = small_n + 1; i < big_n + 1; ++i) {
            result_coefficients[i] = polynomial.coefficients[i];
        }
    }

    return { big_n, result_coefficients };

}

Polynomial Polynomial::operator-() {
    double* result_coefficients;
    result_coefficients = new double[n];
    for (int i = 0; i < n + 1; ++i) {
        result_coefficients[i] = -coefficients[i];
    }
    return { n, result_coefficients };
}

Polynomial Polynomial::operator-(Polynomial polynomial) {

    int big_n = std::max(n, polynomial.n);
    int small_n = std::min(n, polynomial.n);

    double* result_coefficients;
    result_coefficients = new double[big_n];

    for (int i = 0; i < small_n + 1; ++i) {
        result_coefficients[i] = coefficients[i] - polynomial.coefficients[i];
    }

    if (n > small_n) {
        for (int i = small_n + 1; i < big_n + 1; ++i) {
            result_coefficients[i] = coefficients[i];
        }
    }
    else {
        for (int i = small_n + 1; i < big_n + 1; ++i) {
            result_coefficients[i] = -polynomial.coefficients[i];
        }
    }

    return { big_n, result_coefficients };

}

Polynomial Polynomial::operator*(double a) {
    double* result_coefficients;
    result_coefficients = new double[n];
    for (int i = 0; i < n + 1; ++i) {
        result_coefficients[i] = a * coefficients[i];
    }
    return { n, result_coefficients };
}

Polynomial Polynomial::operator*(Polynomial polynomial) {

    auto* result_coefficients = new double[n + polynomial.n + 2];

    for (int i = 0; i < n + polynomial.n + 2; ++i) {
        result_coefficients[i] = 0.;
    }

    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < polynomial.n + 1; ++j) {
            result_coefficients[i + j] += coefficients[i] * polynomial.coefficients[j];
        }
    }

    return { n + polynomial.n, result_coefficients };

}

void Polynomial::present_coefficients() {
    for (int i = 0; i < n + 1; ++i) {
        std::cout << coefficients[i] << " ";
    }
    std::cout << std::endl;
}

Polynomial Polynomial::legendre(int n_) {

    auto* result_coefficients = new double[n_ + 1];

    for (int i = 0; i < n_ + 2; ++i) {
        result_coefficients[i] = 0.;
    }

    int max = 0;

    for (int m = 0; m < div(n_, 2).quot + 1; ++m) {
        result_coefficients[n_ - 2 * m] = pow(-1., m) * factorial(2 * n_ - 2 * m) / (pow(2, n_) * factorial(m) * factorial(n_ - m) * factorial(n_ - 2 * m));
    }

    return { n_, result_coefficients };
}

Polynomial* Polynomial::lagrange(const double* x, int n_) {

    auto* lagrange_polynomials = new Polynomial[n_];
    double denominator;
    double monomial_coefficients[2];
    double unit_coefficient[1] = { 1 };

    for (int i = 0; i < n_; ++i) {
        lagrange_polynomials[i] = Polynomial(0, unit_coefficient);
        for (int j = 0; j < n_; ++j) {
            if (j != i) {
                denominator = x[i] - x[j];
                monomial_coefficients[0] = -x[j] / denominator;
                monomial_coefficients[1] = 1. / denominator;
                lagrange_polynomials[i] = lagrange_polynomials[i] * Polynomial(1, monomial_coefficients);
            }
        }
    }
    return lagrange_polynomials;
}
