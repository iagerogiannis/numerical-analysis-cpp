#include "integration.h"

double integration::trapezoid(double(*f)(double), double x0, double xn, int n)
{
	double I = 0.;
	double dx = (xn - x0) / n;

	for (int i = 0; i < n; ++i) {
		double xa = x0 + i * dx;
		double xb = x0 + (i + 1) * dx;
		double dI = (xb - xa) * (f(xa) + f(xb)) / 2;
		I += dI;
	}
	return I;
}

double integration::simpson1_3(double(*f)(double), double x0, double xn, int n)
{
	double I = 0.;
	double dx = (xn - x0) / n;

	for (int i = 0; i < n; ++i) {
		double xa = x0 + i * dx;
		double xb = x0 + (i + 1) * dx;
		double dI = (xb - xa) * (f(xa) + 4 * f(0.5 * (xa + xb)) + f(xb)) / 6;
		I += dI;
	}
	return I;
}

double integration::simpson3_8(double(*f)(double), double x0, double xn, int n)
{
	double I = 0.;
	double dx = (xn - x0) / n;

	for (int i = 0; i < n; ++i) {
		double xa = x0 + i * dx;
		double xb = x0 + (i + 1) * dx;
		double dI = (xb - xa) * (f(xa) + 3 * f((2 * xa + xb) / 3) +
			3 * f((xa + 2 * xb) / 3) + f(xb)) / 8;
		I += dI;
	}
	return I;
}

double integration::romberg(double(*f)(double), double x0, double xn, int level)
{
	double dx = (xn - x0) / pow(2, level - 1);

	auto* xa = new double[level];

	auto** Irom = new double*[level];
	for (int i = 0; i < level; ++i) {
		Irom[i] = new double[level - i];
		for (int j = 0; j < level - i; ++j) {
			Irom[i][j] = 0.;
		}
	}

	for (int i=1; i < pow(2, level-1)+1; ++i) {
		double xi = x0 + i * dx;
		for (int j = 0; j < level; ++j) {
			if (i % (int)pow(2, level - (j + 1)) == 0) {
				Irom[j][0] += trapezoid(f, xa[j], xi, 1);
				xa[j] = xi;
			}
		}
	}

	for (int j = 1; j < level; ++j) {
		for (int i = 0; i < level - j; ++i) {
			Irom[i][j] = (pow(4, j) * Irom[i + 1][j - 1] - Irom[i][j - 1]) / (pow(4, j) - 1);
		}
	}

	double I = Irom[0][level - 1];

	for (int i = 0; i < level; ++i) {
		delete[] Irom[i];
	}

	delete[] Irom;
	delete[] xa;

	return I;
}

double integration::gauss_legendre(double(*f)(double), double x0, double xn, int n)
{
	return 0.0;
}
