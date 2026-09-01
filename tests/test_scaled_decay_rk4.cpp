#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/rk4.hpp"
#include "odelib/steppers/euler.hpp"
#include "test_utils.hpp"


double k = 0.5;
auto scaled_decay = [](double t, double y) { return k * -y; };

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.05;

    const double exact = y0 * std::exp(-t1 * k);
    odelib::RK4Stepper<double> stepper;

    TestUtils::printHeader("Accuracy Test");
    auto ys1 = odelib::integrate(stepper, scaled_decay, y0, t0, t1, h);
    double error1 = std::abs(ys1.back() - exact);

    TestUtils::printRow("RK4 value (h):", ys1.back());
    TestUtils::printRow("Exact value:", exact);
    TestUtils::printRow("Error (h):", error1);

    assert(error1 < 0.01 && "RK4 solution too far from exact decay solution");
    TestUtils::printPassed("Accuracy Test");

    std::cout << "\nAll test_scaled_decay_rk4 checks passed.\n";

    return 0;
}
