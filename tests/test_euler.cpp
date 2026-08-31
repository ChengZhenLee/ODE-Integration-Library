#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
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
    const double h = 0.001;

    const double exact = y0 * std::exp(-t1);
    odelib::EulerStepper<double> stepper;

    printHeader("Accuracy Test");
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);
    double error = std::abs(ys.back() - exact);

    printRow("Euler value:", ys.back());
    printRow("Exact value:", exact);
    printRow("Error:", error);

    assert(error < 0.01 && "Euler solution too far from exact decay solution");
    printPassed("Accuracy Test");

    std::cout << "\nAll test_euler checks passed.\n";
    return 0;
}
