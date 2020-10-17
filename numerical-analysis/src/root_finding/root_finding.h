#pragma once

#include <functional>

namespace root_finding {

    double secant(double (*f)(double), double x0, double x1, double error = 1e-14);

    double newton_raphson(std::function<double(double)> f, std::function<double(double)> df_dx,
        double x0, double error = 1e-14);

}