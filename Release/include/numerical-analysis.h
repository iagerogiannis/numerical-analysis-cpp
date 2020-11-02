#pragma once

#include <functional>

#include <cpp_extended.h>
#include <Polynomial.h>


namespace root_finding {

    double secant(std::function<double(double)> f, double x0, double x1, double error = 1e-14);

    double newton_raphson(std::function<double(double)> f, std::function<double(double)> df_dx,
        double x0, double error = 1e-14);

    double* newton_raphson_multiple_roots(std::function<double(double)> f, std::function<double(double)> df_dx, int n, double error = 1e-14);

}

namespace integration {

    double trapezoid(std::function<double(double)> f, double x0, double xn, int n);

    double simpson1_3(std::function<double(double)> f, double x0, double xn, int n);

    double simpson3_8(std::function<double(double)> f, double x0, double xn, int n);

    double romberg(std::function<double(double)> f, double x0, double xn, int level);

    double gauss_legendre(std::function<double(double)> f, double x0, double xn, int n);

}

namespace splines {
	
	struct polynomials {
    Polynomial x;
    Polynomial y;
    Polynomial dx_dt;
    Polynomial dy_dt;
    Polynomial* p_i;
	};
	
	struct complex_parameter {
    int s;
    double t;
	};

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

        double x_t(double t);

        double y_t(double t);

        double y_x(double x);

        double dx_dt(double t);

        double dy_dt(double t);

        double dy_dx(double x);

        double p_i(int i, double t) const;

        void presentControlPoints();

        void presentMatrix();

        void presentCoefficients();

    };

    class CompositeQuadraticBezier {
		
    private:
	
        int N{};
        int M{};
        double** ControlPoints{};
        double** EnrichedControlPoints{};
        splines::Bezier* Sectors{};

        void allocateEnrichedCP();
        void allocateSectors();
        double** sectorCP(int);
        static void clearSectorCP(double**);
        void enrichControlPoints();
        void calculateSectors();

    public:
	
        CompositeQuadraticBezier();
        CompositeQuadraticBezier(int N, double** ControlPoints);
        ~CompositeQuadraticBezier();

        void readControlPoints(std::string filename);

        complex_parameter translate_t(double) const;
        double x_t(double);
        double y_t(double);
        double y_x(double);
        double dx_dt(double);
        double dy_dt(double);
    };

}
