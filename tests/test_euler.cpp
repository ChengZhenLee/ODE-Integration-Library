#include <cassert>
#include <cmath>
#include <iostream>

#include "odelib/integrate.hpp"
#include "odelib/steppers/euler.hpp"

double decay(double /*t*/, double y) { return -y; }

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.001;

    odelib::EulerStepper<double> stepper;
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);

    double exact = y0 * std::exp(-t1);
    double error = std::abs(ys.back() - exact);

    std::cout << "Final Euler value: " << ys.back() << "\n";
    std::cout << "Exact value:       " << exact << "\n";
    std::cout << "Error:             " << error << "\n";

    assert(error < 0.01 && "Euler solution too far from exact decay solution");

    std::cout << "test_euler PASSED\n";
    return 0;
}
