#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/adaptive_rk4.hpp"
#include "test_utils.hpp"


double k = 0.5;
auto scaled_decay = [](double t, double y) { return -k * y; };

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.05;

    const double exact = y0 * std::exp(-t1 * k);
    odelib::AdaptiveRK4Stepper<double, double> stepper;

    TestUtils::printHeader("Accuracy Test");
    auto result = odelib::integrate(stepper, scaled_decay, y0, t0, t1, h);
    double error1 = std::abs(result.states.back() - exact);

    TestUtils::printRow("Adaptive RK4 value:", result.states.back());
    TestUtils::printRow("Exact value:", exact);
    TestUtils::printRow("Error:", error1);
    TestUtils::printRow("Final time:", result.times.back());

    assert(error1 < 0.01 && "Adaptive RK4 solution too far from exact decay solution");
    assert(std::abs(result.times.back() - t1) < 1e-9 && "Adaptive RK4 did not integrate to t1");
    TestUtils::printPassed("Accuracy Test");

    std::cout << "\nAll test_scaled_decay_adaptive_rk4 checks passed.\n";

    return 0;
}
