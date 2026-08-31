#include <cassert>
#include <cmath>
#include <iostream>

#include "odelib/integrate.hpp"
#include "odelib/steppers/rk4.hpp"

double decay(double t, double y) { return -y; }

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.001;

    odelib::RK4Stepper<double> stepper;
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);

    double exact = y0 * std::exp(-t1);
    double error = std::abs(ys.back() - exact);

    std::cout << "Final RK4 value: " << ys.back() << "\n";
    std::cout << "Exact value:       " << exact << "\n";
    std::cout << "Error:             " << error << "\n";

    assert(error < 0.01 && "RK4 solution too far from exact decay solution");

    std::cout << "test_rk4 PASSED\n";
    return 0;
}
