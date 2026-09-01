#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/euler.hpp"
#include "test_utils.hpp"


double decay(double t, double y) { return -y; }

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.001;

    const double exact = y0 * std::exp(-t1);
    odelib::EulerStepper<double> stepper;

    TestUtils::printHeader("Accuracy Test");
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);
    double error = std::abs(ys.back() - exact);

    TestUtils::printRow("Euler value:", ys.back());
    TestUtils::printRow("Exact value:", exact);
    TestUtils::printRow("Error:", error);

    assert(error < 0.01 && "Euler solution too far from exact decay solution");
    TestUtils::printPassed("Accuracy Test");

    std::cout << "\nAll test_euler checks passed.\n";
    return 0;
}
