#include "root_finding.h"

double root_finding::secant(double (*f)(double), double x0, double x1, double error) {

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
