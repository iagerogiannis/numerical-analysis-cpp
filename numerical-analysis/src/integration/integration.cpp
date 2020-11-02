#include "integration.h"

#include <iostream>
#include <functional>
#include "Polynomial.h"
#include "..\root_finding\root_finding.h"

using namespace std::placeholders;

double integration::trapezoid(std::function<double(double)> f, double x0, double xn, int n)
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

double integration::simpson1_3(std::function<double(double)> f, double x0, double xn, int n)
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

double integration::simpson3_8(std::function<double(double)> f, double x0, double xn, int n)
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

double integration::romberg(std::function<double(double)> f, double x0, double xn, int level)
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

double integration::gauss_legendre(std::function<double(double)> f, double x0, double xn, int n)
{
	Polynomial p = Polynomial::legendre(n + 1);
	
	auto p_x = std::bind(&Polynomial::value, p, _1);
	auto dp_dx = std::bind(&Polynomial::value, p.derivative(), _1);

	double* x = root_finding::newton_raphson_multiple_roots(p_x, dp_dx, n + 1);

	Polynomial* L = Polynomial::lagrange(x, n + 1);

	auto w = new double[n + 1];
	for (int i = 0; i < n + 1; ++i) {
		w[i] = L[i].definite_integral(-1., 1.);
	}

	double I = 0.;
	for (int i = 0; i < n + 1; ++i) {
		I += w[i] * f(0.5 * (x[i] * (xn - x0) + xn + x0));
	}

	I *= 0.5 * (xn - x0);

	delete[] w;

	return I;
}
