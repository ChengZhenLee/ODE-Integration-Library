#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/rk4.hpp"
#include "odelib/steppers/euler.hpp"

double decay(double t, double y) { return -y; }

namespace {

void printHeader(const std::string& title) {
    std::cout << "\n-- " << title << " " << std::string(std::max<int>(0, 40 - static_cast<int>(title.size())), '-') << "\n";
}

void printRow(const std::string& label, double value) {
    std::cout << "  " << std::left << std::setw(22) << label
               << std::right << std::scientific << std::setprecision(6) << value << "\n";
}

void printPassed(const std::string& name) {
    std::cout << "[PASSED] " << name << "\n";
}

} 

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.05;

    const double exact = y0 * std::exp(-t1);
    odelib::RK4Stepper<double> stepper;

    printHeader("Accuracy Test");
    auto ys1 = odelib::integrate(stepper, decay, y0, t0, t1, h);
    double error1 = std::abs(ys1.back() - exact);

    printRow("RK4 value (h):", ys1.back());
    printRow("Exact value:", exact);
    printRow("Error (h):", error1);

    assert(error1 < 0.01 && "RK4 solution too far from exact decay solution");
    printPassed("Accuracy Test");

    printHeader("Convergence Test");
    auto ys2 = odelib::integrate(stepper, decay, y0, t0, t1, h / 2);
    double error2 = std::abs(ys2.back() - exact);
    double errorRatio = error1 / error2;

    printRow("RK4 value (h/2):", ys2.back());
    printRow("Error (h/2):", error2);
    printRow("Error ratio:", errorRatio);

    assert(errorRatio > 14 && errorRatio < 18 && "RK4 error ratio not close to expected 16x");
    printPassed("Convergence Test");

    printHeader("Euler Comparison Test");
    odelib::EulerStepper<double> eulerStepper;
    auto ys3 = odelib::integrate(eulerStepper, decay, y0, t0, t1, h);
    double error3 = std::abs(ys3.back() - exact);

    printRow("Euler error:", error3);
    printRow("RK4 error:", error1);

    assert(error1 < error3 / 100 && "RK4 error should be dramatically smaller than Euler's error");
    printPassed("Euler Comparison Test");

    std::cout << "\nAll test_rk4 checks passed.\n";
    return 0;
}
