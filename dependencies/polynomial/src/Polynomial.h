#pragma once

class Polynomial {

private:
    double* coefficients{};
    int n{};

public:

    Polynomial();
    Polynomial(int n, double* coefficients);
    Polynomial(const Polynomial&);
    Polynomial(Polynomial&&) noexcept;
    Polynomial& operator = (Polynomial&&) noexcept;

    ~Polynomial();

    void setCoefficients(int, double*);

    double value(double x);
    double definite_integral(double a, double b);

    Polynomial derivative();
    Polynomial indefinite_integral();

    Polynomial operator + (Polynomial);

    Polynomial operator - ();
    Polynomial operator - (Polynomial);

    Polynomial operator * (double);
    Polynomial operator * (Polynomial);

    static Polynomial legendre(int);

    static Polynomial* lagrange(int n_, const double* x);

    void present_coefficients();

};
