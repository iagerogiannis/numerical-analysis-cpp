#pragma once

#include <functional>

namespace integration {

    double trapezoid(std::function<double(double)> f, double x0, double xn, int n);

    double simpson1_3(std::function<double(double)> f, double x0, double xn, int n);

    double simpson3_8(std::function<double(double)> f, double x0, double xn, int n);

    double romberg(std::function<double(double)> f, double x0, double xn, int level);

    double gauss_legendre(std::function<double(double)> f, double x0, double xn, int n);

}
