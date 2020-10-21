#pragma once

#include <functional>

namespace integration {

    double trapezoid(double (*f)(double), double x0, double xn, int n);

    double simpson1_3(double (*f)(double), double x0, double xn, int n);

    double simpson3_8(double (*f)(double), double x0, double xn, int n);

    double romberg(double (*f)(double), double x0, double xn, int level);

    double gauss_legendre(double (*f)(double), double x0, double xn, int n);

}
