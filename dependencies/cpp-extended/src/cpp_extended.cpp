#include "cpp_extended.h"

#include <iostream>

double* arrays::copyDynamicArray(const double* initial_array, int N) {

    auto* destination = new double[N];

    for (int i = 0; i < N; ++i) {
        destination[i] = initial_array[i];
    }

    return destination;
}

double** arrays::copyDynamicArray(double** initial_array, int N, int M) {

    auto** destination = new double* [N];
    for (int i = 0; i < N; ++i) {
        destination[i] = new double[M];
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            destination[i][j] = initial_array[i][j];
        }
    }

    return destination;
}

void arrays::quickSort(double* A, int in, bool ascending, int i0)
{
    if (i0 < in)
    {
        int p = partition(A, i0, in, ascending);
        quickSort(A, p, ascending, i0);
        quickSort(A, in, ascending, p + 1);
    }
}

int arrays::partition(double* A, int i0, int in, bool ascending)
{
    int p = i0;

    if (ascending) {
        for (int j = i0 + 1; j < in; ++j)
        {
            if (A[j] <= A[i0])
            {
                p++;
                std::swap(A[p], A[j]);
            }
        }
    }
    else {
        for (int j = i0 + 1; j < in; ++j)
        {
            if (A[j] >= A[i0])
            {
                p++;
                std::swap(A[p], A[j]);
            }
        }
    }

    std::swap(A[p], A[i0]);
    return p;
}

long long int math::factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
