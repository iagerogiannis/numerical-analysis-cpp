#include "Bezier.h"

#include <fstream>
#include <functional>
#include <cmath>

#include "cpp_extended.h"
#include "Polynomial.h"
#include "..\..\root_finding\root_finding.h"

using namespace std::placeholders;

splines::Bezier::Bezier() = default;

splines::Bezier::Bezier(int N, double **cp)
: N(N)
{
    ControlPoints = arrays::copyDynamicArray(cp, 2, N + 1);

    allocateMatrix();
    allocateCoefficients();
    allocatePolynomials();

    calculateMatrix();
    calculateCoefficients();

    Polynomial x_polynomial(N, Coefficients[0]);
    Polynomial y_polynomial(N, Coefficients[1]);

    Polynomials.x = std::move(x_polynomial);
    Polynomials.y = std::move(y_polynomial);

    Polynomials.dx_dt = Polynomials.x.derivative();
    Polynomials.dy_dt = Polynomials.y.derivative();

    for (int i=0; i < N + 1; ++i) {
        Polynomial p_i(N, Matrix[i]);
        Polynomials.p_i[i] = std::move(p_i);
    }
}

splines::Bezier::Bezier(const Bezier &bezier) {

    N = bezier.N;
    ControlPoints = arrays::copyDynamicArray(bezier.ControlPoints, 2, N + 1);
    Matrix = arrays::copyDynamicArray(bezier.Matrix, N + 1, N + 1);
    Coefficients = arrays::copyDynamicArray(bezier.ControlPoints, 2, N + 1);

    allocatePolynomials();
    copyPolynomials(bezier);
}

splines::Bezier::Bezier(Bezier &&bezier) noexcept {

    N = bezier.N;
    ControlPoints = bezier.ControlPoints;
    Matrix = bezier.Matrix;
    Coefficients = bezier.Coefficients;

    Polynomials.x = std::move(bezier.Polynomials.x);
    Polynomials.y = std::move(bezier.Polynomials.y);
    Polynomials.dx_dt = std::move(bezier.Polynomials.dx_dt);
    Polynomials.dy_dt = std::move(bezier.Polynomials.dy_dt);
    Polynomials.p_i = bezier.Polynomials.p_i;

    bezier.ControlPoints = nullptr;
    bezier.Matrix = nullptr;
    bezier.Coefficients = nullptr;
    bezier.Polynomials.p_i = nullptr;
}

splines::Bezier &splines::Bezier::operator=(Bezier &&bezier)  noexcept {
    if (this != &bezier) {
        if (ControlPoints != nullptr) {
            for (int i=0; i < 2; ++i) {
                delete[] ControlPoints[i];
            }
        }
        if (Matrix != nullptr) {
            for (int i=0; i < N + 1; ++i) {
                delete[] Matrix[i];
            }
        }
        if (Coefficients != nullptr) {
            for (int i=0; i < 2; ++i) {
                delete[] Coefficients[i];
            }
        }

        delete[] ControlPoints;
        delete[] Matrix;
        delete[] Coefficients;
        delete[] Polynomials.p_i;

        N = bezier.N;
        ControlPoints = bezier.ControlPoints;
        Matrix = bezier.Matrix;
        Coefficients = bezier.Coefficients;

        Polynomials.x = std::move(bezier.Polynomials.x);
        Polynomials.y = std::move(bezier.Polynomials.y);
        Polynomials.dx_dt = std::move(bezier.Polynomials.dx_dt);
        Polynomials.dy_dt = std::move(bezier.Polynomials.dy_dt);
        Polynomials.p_i = bezier.Polynomials.p_i;

        bezier.ControlPoints = nullptr;
        bezier.Matrix = nullptr;
        bezier.Coefficients = nullptr;
        bezier.Polynomials.p_i = nullptr;
    }
    return *this;
}

splines::Bezier::~Bezier() {

    if (ControlPoints != nullptr) {
        for (int i=0; i < 2; ++i) {
            delete[] ControlPoints[i];
        }
    }
    if (Matrix != nullptr) {
        for (int i=0; i < N + 1; ++i) {
            delete[] Matrix[i];
        }
    }
    if (Coefficients != nullptr) {
        for (int i=0; i < 2; ++i) {
            delete[] Coefficients[i];
        }
    }

    delete[] ControlPoints;
    delete[] Matrix;
    delete[] Coefficients;
    delete[] Polynomials.p_i;

    ControlPoints = nullptr;
    Matrix = nullptr;
    Coefficients = nullptr;
    Polynomials.p_i = nullptr;
}

void splines::Bezier::allocatePolynomials() {
    Polynomials.p_i = new Polynomial[N + 1];
}

void splines::Bezier::allocateMatrix() {
    Matrix = new double*[N + 1];
    for (int i=0; i < N + 1; ++i) {
        Matrix[i] = new double[N + 1];
    }
}

void splines::Bezier::allocateCoefficients() {
    Coefficients = new double*[2];
    for (int i=0; i < 2; ++i) {
        Coefficients[i] = new double[N + 1];
        for (int j=0; j < N + 1; ++j) {
            Coefficients[i][j] = 0.;
        }
    }
}

void splines::Bezier::copyPolynomials(const Bezier &bezier) {
    {
        Polynomial p = bezier.Polynomials.x;
        Polynomials.x = std::move(p);
    }
    {
        Polynomial p = bezier.Polynomials.y;
        Polynomials.y = std::move(p);
    }
    {
        Polynomial p = bezier.Polynomials.dx_dt;
        Polynomials.dx_dt = std::move(p);
    }
    {
        Polynomial p = bezier.Polynomials.dy_dt;
        Polynomials.dy_dt = std::move(p);
    }

    for (int i = 0; i < N + 1; ++i) {
        Polynomial p = bezier.Polynomials.p_i[i];
        Polynomials.p_i[i] = std::move(p);
    }
}

void splines::Bezier::readControlPoints(std::string filename) {

    ControlPoints = new double* [2];

    std::ifstream cp_file;
    std::string line;

    cp_file.open(filename);
    
    int i = 0;
    while (std::getline(cp_file, line)) {
        if (i == 0) {
            N = std::stoi(line) - 1;
            ControlPoints[0] = new double[N + 1];
            ControlPoints[1] = new double[N + 1];
        } 
        else {
            ControlPoints[div((i - 1), (N + 1)).quot][(i - 1) % (N + 1)] = std::stod(line);
        }
        ++i;
    }

    cp_file.close();

    allocateMatrix();
    allocateCoefficients();
    allocatePolynomials();

    calculateMatrix();
    calculateCoefficients();

    Polynomial x_polynomial(N, Coefficients[0]);
    Polynomial y_polynomial(N, Coefficients[1]);

    Polynomials.x = std::move(x_polynomial);
    Polynomials.y = std::move(y_polynomial);

    Polynomials.dx_dt = Polynomials.x.derivative();
    Polynomials.dy_dt = Polynomials.y.derivative();

    for (int i = 0; i < N + 1; ++i) {
        Polynomial p_i(N, Matrix[i]);
        Polynomials.p_i[i] = std::move(p_i);
    }
}

void splines::Bezier::calculateMatrix() {
    for (int i=0; i < N + 1; ++i) {
        for (int j=0; j < N + 1; ++j) {
            if (j<i) {
                Matrix[i][j] = 0.;
            } else {
                Matrix[i][j] = (double)(pow(-1, j-i) * math::factorial(N) / (math::factorial(i) * math::factorial(j-i) * math::factorial(N - j)));
            }
        }
    }

}

void splines::Bezier::calculateCoefficients() {
    for (int k=0; k < 2; ++k) {
        for (int i=0; i < N + 1; ++i) {
            for (int j=0; j < N + 1; ++j) {
                Coefficients[k][i] += Matrix[j][i] * ControlPoints[k][j];
            }
        }
    }
}

double splines::Bezier::x_t(double t) {
    return Polynomials.x.value(t);
}

double splines::Bezier::y_t(double t) {
    return Polynomials.y.value(t);
}

double splines::Bezier::y_x(double x) {
    auto x_ = std::bind(&splines::Bezier::x_t, this, _1);
    auto dx_dt_ = std::bind(&splines::Bezier::dx_dt, this, _1);
    double t0 = root_finding::newton_raphson([&x_, &x](double t) {return x_(t) - x; }, dx_dt_, 0.5);
    return y_t(t0);
}

double splines::Bezier::dx_dt(double t) {
    return Polynomials.dx_dt.value(t);
}

double splines::Bezier::dy_dt(double t) {
    return Polynomials.dy_dt.value(t);
}

double splines::Bezier::dy_dx(double x) {
    auto x_ = std::bind(&splines::Bezier::x_t, this, _1);
    auto dx_dt_ = std::bind(&splines::Bezier::dx_dt, this, _1);
    double t0 = root_finding::newton_raphson([&x_, &x](double t) {return x_(t) - x; }, dx_dt_, 0.5);
    return dy_dt(t0) / dx_dt(t0);
}

double splines::Bezier::p_i(int i, double t) const {
    return Polynomials.p_i[i].value(t);
}

void splines::Bezier::presentControlPoints(){
    for (int i=0; i < N + 1; ++i) {
        std::cout << ControlPoints[0][i] << " " << ControlPoints[1][i] << std::endl;
    }
}

void splines::Bezier::presentMatrix(){
    for (int i=0; i < N + 1; ++i) {
        for (int j=0; j < N + 1; ++j) {
            std::cout << Matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void splines::Bezier::presentCoefficients() {
    for (int i=0; i < 2; ++i) {
        for (int j=0; j < N + 1; ++j) {
            std::cout << Coefficients[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
