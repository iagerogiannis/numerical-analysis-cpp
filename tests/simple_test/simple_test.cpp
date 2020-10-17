#include <iostream>
#include <functional>
#include <chrono>

#include "numerical-analysis.h"

using namespace std::chrono;

double f(double x) {
    return pow(x-1, 2) - 3 * x;
}

double df_dx(double x) {
    return 2 * x - 5;
}

int main() {

    std::cout << root_finding::newton_raphson(f, df_dx, 0) << std::endl;

    int n = 4;

    auto** cp = new double* [2];
    for (int i = 0; i < 2; ++i) {
        cp[i] = new double[n + 1];
    }

    cp[0][0] = 1.;
    cp[0][1] = 2.;
    cp[0][2] = 3.;
    cp[0][3] = 4.;
    cp[0][4] = 5.;
    cp[1][0] = 3.;
    cp[1][1] = 4.5;
    cp[1][2] = 5.;
    cp[1][3] = 1.;
    cp[1][4] = 1.;

    splines::CompositeQuadraticBezier comp_bezier(n + 1, cp);

    int points = 1000;

    auto start = high_resolution_clock::now();

    for (int i = 0; i < points + 1; ++i) {
        double x = (1. / points) * (double)i * 4. + 1.;
        double y = comp_bezier.y_x(x);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "time taken by function: " << 1e-3 * duration.count() << " milliseconds" << std::endl;

    std::cin.get();

}
