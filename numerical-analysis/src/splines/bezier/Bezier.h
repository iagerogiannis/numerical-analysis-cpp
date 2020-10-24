#pragma once

#include <iostream>
#include <string>
#include "Polynomial.h"

struct polynomials {
    Polynomial x;
    Polynomial y;
    Polynomial dx_dt;
    Polynomial dy_dt;
    Polynomial* p_i;
};

namespace splines {

    class Bezier {

    private:
        int N{};
        double** ControlPoints{};
        double** Matrix{};
        double** Coefficients{};

        polynomials Polynomials{};

        void allocateMatrix();
        void allocateCoefficients();
        void allocatePolynomials();
        void copyPolynomials(const Bezier&);

        void calculateMatrix();
        void calculateCoefficients();

    public:

        Bezier();
        Bezier(int N, double** cp);
        Bezier(const Bezier&);
        Bezier(Bezier&&) noexcept;
        Bezier& operator = (Bezier&&) noexcept;


        ~Bezier();

        void readControlPoints(std::string filename);

        double y_x(double x);

        double x_t(double t);

        double y_t(double t);

        double dx_dt(double t);

        double dy_dt(double t);

        double p_i(int i, double t) const;

        void presentControlPoints();

        void presentMatrix();

        void presentCoefficients();

    };
}
