#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/rkdp.hpp"
#include "Eigen/Dense"
#include "test_utils.hpp"


auto harmonic_oscillator = [](double t, Eigen::Vector2d y) { return Eigen::Vector2d(y(1), -y(0)); };

int main() {
    const Eigen::Vector2d y0(1.0, 2.0);
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.05;

    const Eigen::Vector2d exact = Eigen::Vector2d(y0(0) * cos(t1) + y0(1) * sin(t1), -y0(0) * sin(t1) + y0(1) * cos(t1));
    odelib::RKDPStepper<Eigen::Vector2d, double> stepper;

    TestUtils::printHeader("Accuracy Test");
    auto ys1 = odelib::integrate(stepper, harmonic_oscillator, y0, t0, t1, h);
    Eigen::Vector2d error1 = (ys1.back() - exact).cwiseAbs();

    TestUtils::printRow("RKDP value (h):", ys1.back().transpose());
    TestUtils::printRow("Exact value:", exact.transpose());
    TestUtils::printRow("Error (h):", error1.transpose());

    assert(error1(0) < 0.01 && error1(1) < 0.01 && "RKDP solution too far from exact harmonic oscillator solution");
    TestUtils::printPassed("Accuracy Test");

    std::cout << "\nAll test_harmonic_oscillator_rkdp checks passed.\n";

    return 0;
}
