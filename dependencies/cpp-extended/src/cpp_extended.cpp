#include "cpp_extended.h"

double* copyDynamicArray(const double* original, int N) {

    auto* destination = new double[N];

    for (int i = 0; i < N; ++i) {
        destination[i] = original[i];
    }

    return destination;
}

double** copyDynamicArray(double** original, int N, int M) {

    auto** destination = new double* [N];
    for (int i = 0; i < N; ++i) {
        destination[i] = new double[M];
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            destination[i][j] = original[i][j];
        }
    }

    return destination;
}
