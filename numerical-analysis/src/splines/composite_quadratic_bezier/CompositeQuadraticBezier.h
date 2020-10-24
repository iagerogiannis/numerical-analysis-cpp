#pragma once

#include <iostream>
#include <string>
#include "..\bezier\Bezier.h"

struct complex_parameter {
    int s;
    double t;
};

namespace splines {
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
