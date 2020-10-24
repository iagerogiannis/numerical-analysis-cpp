#include <iostream>
#include <functional>
#include <chrono>
#include <string>

#include "numerical-analysis.h"

using namespace std::chrono;

double f(double x) {
    return pow(x-1, 2) - 3 * x;
}

double df_dx(double x) {
    return 2 * x - 5;
}

int main() {

    std::cout << root_finding::secant(f, 0., 1.) << std::endl;
    std::cout << root_finding::newton_raphson(f, df_dx, 0) << std::endl;

    splines::Bezier curve;
    curve.readControlPoints("control_points.dat");

    int points = 1000;

    for (int i = 0; i < points + 1; ++i) {
        double t = (double)i / points;
        double x = curve.x_t(t);
        double y = curve.y_t(t);
        std::cout << x << " " << y << "\n";
    }

    double x0 = 0.;
    double xn = 10.;
    int n_ = 1000000;
    int rom_level = 10;
    int gl_level = 8;

    auto start = high_resolution_clock::now();

    std::cout << integration::trapezoid(f, x0, xn, n_) << std::endl;
    std::cout << integration::simpson1_3(f, x0, xn, n_) << std::endl;
    std::cout << integration::simpson3_8(f, x0, xn, n_) << std::endl;
    std::cout << integration::romberg(f, x0, xn, rom_level) << std::endl;
    std::cout << integration::gauss_legendre(f, x0, xn, gl_level) << std::endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    std::cout << "Time needed by function: " << 1e-3 * duration.count() << " milliseconds" << std::endl;

    std::cout << std::endl << "Press ENTER to Continue..." << std::endl;
    std::cin.get();

}
