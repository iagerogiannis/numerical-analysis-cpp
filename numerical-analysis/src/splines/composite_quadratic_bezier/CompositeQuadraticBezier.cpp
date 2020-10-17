#include "CompositeQuadraticBezier.h"

#include <iostream>
#include "cpp_extended.h"

splines::CompositeQuadraticBezier::CompositeQuadraticBezier() = default;

splines::CompositeQuadraticBezier::CompositeQuadraticBezier(int N, double** ControlPoints_)
    : N(N) {
    M = 2 * N - 3;
    allocateEnrichedCP();
    allocateSectors();
    ControlPoints = copyDynamicArray(ControlPoints_, 2, N);
    enrichControlPoints();
    calculateSectors();
}

splines::CompositeQuadraticBezier::~CompositeQuadraticBezier() {
    if (ControlPoints != nullptr) {
        for (int i = 0; i < 2; ++i) {
            delete[] ControlPoints[i];
        }
    }
    if (EnrichedControlPoints != nullptr) {
        for (int i = 0; i < 2; ++i) {
            delete[] EnrichedControlPoints[i];
        }
    }

    delete[] ControlPoints;
    delete[] EnrichedControlPoints;
    delete[] Sectors;

    ControlPoints = nullptr;
    EnrichedControlPoints = nullptr;
    Sectors = nullptr;
}

void splines::CompositeQuadraticBezier::enrichControlPoints() {

    EnrichedControlPoints[0][0] = ControlPoints[0][0];
    EnrichedControlPoints[1][0] = ControlPoints[1][0];

    for (int i = 1; i < M - 1; i += 2) {
        EnrichedControlPoints[0][i] = ControlPoints[0][(i + 1) / 2];
        EnrichedControlPoints[1][i] = ControlPoints[1][(i + 1) / 2];
    }

    for (int i = 2; i < M - 2; i += 2) {
        EnrichedControlPoints[0][i] = 0.5 * (ControlPoints[0][i / 2] + ControlPoints[0][i / 2 + 1]);
        EnrichedControlPoints[1][i] = 0.5 * (ControlPoints[1][i / 2] + ControlPoints[1][i / 2 + 1]);
    }

    EnrichedControlPoints[0][M - 1] = ControlPoints[0][N - 1];
    EnrichedControlPoints[1][M - 1] = ControlPoints[1][N - 1];
}

void splines::CompositeQuadraticBezier::allocateEnrichedCP() {
    EnrichedControlPoints = new double* [2];
    for (int i = 0; i < 2; ++i) {
        EnrichedControlPoints[i] = new double[M];
    }
}

void splines::CompositeQuadraticBezier::allocateSectors() {
    Sectors = new Bezier[N - 2];
}

void splines::CompositeQuadraticBezier::calculateSectors() {
    for (int i = 0; i < N - 2; ++i) {
        double** cp = sectorCP(i);
        Bezier sector(2, cp);
        Sectors[i] = std::move(sector);
        clearSectorCP(cp);
    }
}

double** splines::CompositeQuadraticBezier::sectorCP(int s) {
    auto** cp = new double* [2];
    for (int i = 0; i < 2; ++i) {
        cp[i] = new double[3];
        cp[i][0] = EnrichedControlPoints[i][2 * s];
        cp[i][1] = EnrichedControlPoints[i][2 * s + 1];
        cp[i][2] = EnrichedControlPoints[i][2 * s + 2];
    }
    return cp;
}

void splines::CompositeQuadraticBezier::clearSectorCP(double** cp) {
    delete[] cp[0];
    delete[] cp[1];
}

complex_parameter splines::CompositeQuadraticBezier::translate_t(double t) const {
    complex_parameter translation{};
    if (t != 1.) {
        translation.s = (int)(t * (N - 2));
        translation.t = t * (N - 2) - translation.s;
    }
    else {
        translation.s = N - 3;
        translation.t = 1.;
    }
    return translation;
}

double splines::CompositeQuadraticBezier::y_x(double x) {
    int s = 0;
    while (true) {
        if (Sectors[s].x_t(1.) >= x) {
            break;
        }
        ++s;
    }
    return Sectors[s].y_x(x);
}

double splines::CompositeQuadraticBezier::x_t(double t) {
    complex_parameter ct = translate_t(t);
    return Sectors[ct.s].x_t(ct.t);
}

double splines::CompositeQuadraticBezier::y_t(double t) {
    complex_parameter ct = translate_t(t);
    return Sectors[ct.s].y_t(ct.t);
}

double splines::CompositeQuadraticBezier::dx_dt(double t) {
    complex_parameter ct = translate_t(t);
    return Sectors[ct.s].dx_dt(ct.t);
}

double splines::CompositeQuadraticBezier::dy_dt(double t) {
    complex_parameter ct = translate_t(t);
    return Sectors[ct.s].dy_dt(ct.t);
}
