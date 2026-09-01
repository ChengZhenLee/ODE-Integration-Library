#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/adaptive_rk4.hpp"
#include "Eigen/Dense"
#include "test_utils.hpp"


auto harmonic_oscillator = [](double t, Eigen::Vector2d y) { return Eigen::Vector2d(y(1), -y(0)); };

int main() {
    const Eigen::Vector2d y0(1.0, 2.0);
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.05;

    const Eigen::Vector2d exact = Eigen::Vector2d(y0(0) * cos(t1) + y0(1) * sin(t1), -y0(0) * sin(t1) + y0(1) * cos(t1));
    odelib::AdaptiveRK4Stepper<Eigen::Vector2d, double> stepper;

    TestUtils::printHeader("Accuracy Test");
    auto result = odelib::integrate(stepper, harmonic_oscillator, y0, t0, t1, h);
    Eigen::Vector2d error1 = (result.states.back() - exact).cwiseAbs();

    TestUtils::printRow("Adaptive RK4 value:", result.states.back().transpose());
    TestUtils::printRow("Exact value:", exact.transpose());
    TestUtils::printRow("Error:", error1.transpose());
    TestUtils::printRow("Final time:", result.times.back());

    assert(error1(0) < 0.01 && error1(1) < 0.01 && "Adaptive RK4 solution too far from exact harmonic oscillator solution");
    assert(std::abs(result.times.back() - t1) < 1e-9 && "Adaptive RK4 did not integrate to t1");
    TestUtils::printPassed("Accuracy Test");

    std::cout << "\nAll test_harmonic_oscillator_adaptive_rk4 checks passed.\n";

    return 0;
}
