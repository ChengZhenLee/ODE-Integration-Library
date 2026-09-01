#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/rk4.hpp"
#include "odelib/steppers/euler.hpp"
#include "test_utils.hpp"

double decay(double t, double y) { return -y; }

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.05;

    const double exact = y0 * std::exp(-t1);
    odelib::RK4Stepper<double> stepper;

    TestUtils::printHeader("Accuracy Test");
    auto ys1 = odelib::integrate(stepper, decay, y0, t0, t1, h);
    double error1 = std::abs(ys1.back() - exact);

    TestUtils::printRow("RK4 value (h):", ys1.back());
    TestUtils::printRow("Exact value:", exact);
    TestUtils::printRow("Error (h):", error1);

    assert(error1 < 0.01 && "RK4 solution too far from exact decay solution");
    TestUtils::printPassed("Accuracy Test");

    TestUtils::printHeader("Convergence Test");
    auto ys2 = odelib::integrate(stepper, decay, y0, t0, t1, h / 2);
    double error2 = std::abs(ys2.back() - exact);
    double errorRatio = error1 / error2;

    TestUtils::printRow("RK4 value (h/2):", ys2.back());
    TestUtils::printRow("Error (h/2):", error2);
    TestUtils::printRow("Error ratio:", errorRatio);

    assert(errorRatio > 14 && errorRatio < 18 && "RK4 error ratio not close to expected 16x");
    TestUtils::printPassed("Convergence Test");

    TestUtils::printHeader("Euler Comparison Test");
    odelib::EulerStepper<double> eulerStepper;
    auto ys3 = odelib::integrate(eulerStepper, decay, y0, t0, t1, h);
    double error3 = std::abs(ys3.back() - exact);

    TestUtils::printRow("Euler error:", error3);
    TestUtils::printRow("RK4 error:", error1);

    assert(error1 < error3 / 100 && "RK4 error should be dramatically smaller than Euler's error");
    TestUtils::printPassed("Euler Comparison Test");

    std::cout << "\nAll test_rk4 checks passed.\n";
    return 0;
}
