#include <iostream>
#include <functional>
#include <chrono>
#include <string>

#include "numerical-analysis.h"

using namespace std::placeholders;
using namespace std::chrono;

void measurePerformance(std::function<void(void)> f) {
    auto start = high_resolution_clock::now();

    f();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    std::cout << "Time needed by function: " << 1e-3 * duration.count() << " milliseconds" << std::endl;
}

int main() {

    splines::Bezier curve;
    curve.readControlPoints("control_points.dat");

    int points = 1000;

    for (int i = 0; i < points + 1; ++i) {
        double t = (double)i / points;
        double x = curve.x_t(t);
        double y = curve.y_t(t);
        //std::cout << x << " " << y << "\n";
    }

    auto f = std::bind(&splines::Bezier::y_x, curve, _1);
    auto df_dx = std::bind(&splines::Bezier::dy_dx, curve, _1);

    double x0 = 0.;
    double x1 = 1.;
    double x2 = 3.;
    double xn = 5.;
    int n_ = 1e6;
    int rom_level = 12;
    int gl_level = 8;

    //measurePerformance([&f, &x0, &x1]() {
    //    std::cout << root_finding::secant(f, x0, x1) << std::endl; });
    //measurePerformance([&f, &df_dx, &x2]() {
    //    std::cout << root_finding::newton_raphson(f, df_dx, x2) << std::endl; });

    //measurePerformance([&f, &x0, &xn, &n_]() {
    //    std::cout << integration::trapezoid(f, x0, xn, n_) << std::endl; });
    //measurePerformance([&f, &x0, &xn, &n_]() {
    //    std::cout << integration::simpson1_3(f, x0, xn, n_) << std::endl; });
    //measurePerformance([&f, &x0, &xn, &n_]() {
    //    std::cout << integration::simpson3_8(f, x0, xn, n_) << std::endl; });
    //measurePerformance([&f, &x0, &xn, &rom_level]() {
    //    std::cout << integration::romberg(f, x0, xn, rom_level) << std::endl; });
    //measurePerformance([&f, &x0, &xn, &gl_level](){
    //    std::cout << integration::gauss_legendre(f, x0, xn, gl_level) << std::endl; });
    
    std::cout << std::endl << "Press ENTER to Continue..." << std::endl;
    std::cin.get();
}
