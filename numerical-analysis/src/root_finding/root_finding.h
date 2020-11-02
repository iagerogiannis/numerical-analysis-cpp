#pragma once

#include <functional>

namespace root_finding {

    double secant(std::function<double(double)> f, double x0, double x1, double error = 1e-14);

    double newton_raphson(std::function<double(double)> f, std::function<double(double)> df_dx,
        double x0, double error = 1e-14);

    double* newton_raphson_multiple_roots(std::function<double(double)> f, std::function<double(double)> df_dx, int n, double error = 1e-14);

}
