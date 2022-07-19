#include "ODEIntegrator.hpp"
#include <cmath>
#include <iostream>

struct Complex
{
    double real;
    double imag;
    friend auto operator*(Complex z, double x)
    {
        return Complex{z.real*x, z.imag*x};
    }
    friend auto operator*(double x, Complex z)
    {
        return Complex{z.real*x, z.imag*x};
    }
    friend auto operator+(Complex z1, Complex z2)
    {
        return Complex{z1.real + z2.real, z1.imag + z2.imag};
    }
    friend auto operator/(Complex z, double x)
    {
        return Complex{z.real / x, z.imag / x};
    }
    friend auto operator<<(std::ostream& os, Complex z) -> std::ostream&
    {
        return os << z.real << "+" << z.imag << "i";
    }
};

auto f(double x) -> Complex
{
    return Complex{ cos(x), sin(x) };
}

auto df_dx(const Complex& z) -> Complex
{
    return Complex{ -z.imag, z.real };
}

int main()
{
    auto euler_integrator = ODEIntegrator<Complex, EulerStepper>{};
    auto rk2_integrator   = ODEIntegrator<Complex, RK2Stepper>{};
    auto rk4_integrator   = ODEIntegrator<Complex, RK4Stepper>{};
    auto initial_condition = Complex{1.0, 0.0};
    auto euler_result = euler_integrator.solve(initial_condition, df_dx, 0.0, 10.0);
    auto rk2_result   = rk2_integrator.solve(initial_condition, df_dx, 0.0, 10.0);
    auto rk4_result   = rk4_integrator.solve(initial_condition, df_dx, 0.0, 10.0);
    auto actual       = f(10.0);
    std::cout << "Euler result: " << euler_result << '\n'
              << "RK2 result: "   << rk2_result   << '\n'
              << "RK4 result: "   << rk4_result   << '\n'
              << "Actual value: " << actual       << '\n';
    return 0;
}