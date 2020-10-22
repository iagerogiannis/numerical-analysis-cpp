#include "root_finding.h"

#include <iostream>
#include "cpp_extended.h"

double root_finding::secant(std::function<double(double)> f, double x0, double x1, double error) {

    double fx0 = f(x0);
    double fx1 = f(x1);
    while (std::abs(fx1) > error) {
        double x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0);
        x0 = x1;
        x1 = x2;
        fx0 = f(x0);
        fx1 = f(x2);
    }
    return x1;
}

double root_finding::newton_raphson(std::function<double(double)> f, std::function<double(double)> df_dx,
    double x0, double error) {
    while (std::abs(f(x0)) > error) {
        x0 -= f(x0) / df_dx(x0);
    }
    return x0;
}

double* root_finding::newton_raphson_multiple_roots(std::function<double(double)> f, std::function<double(double)> df_dx, int n, double error)
{
    auto* roots = new double[n];
    for (int i = 0; i < n; ++i) {
        double xk_old = 2.;
        double xk_new;
        while (true) {
            double sigma = 0.;
            for (int j = 0; j < i; ++j) {
                sigma += 1. / (xk_old - roots[j]);
            }
            xk_new = xk_old - f(xk_old) / (df_dx(xk_old) - f(xk_old) * sigma);
            if (abs(xk_new - xk_old) < error) {
                break;
            }
            xk_old = xk_new;
        }
        roots[i] = xk_new;
    }
    quickSort(roots, n);
    return roots;
}
