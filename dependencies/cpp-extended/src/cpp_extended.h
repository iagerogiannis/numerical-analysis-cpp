#pragma once

double* copyDynamicArray(const double* initial_array, int N);
double** copyDynamicArray(double** initial_array, int N, int M);
void quickSort(double* myArray, int N, bool ascending = true, int i0 = 0);
int partition(double* myArray, int i0, int in, bool ascending = true);
