#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/euler.hpp"
#include "odelib/integrate.hpp"

double decay(double t, double y) { return -y; }

int main() {
    const double y0 = 1.0;   // initial condition
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;   // step size

    odelib::EulerStepper<double> stepper;
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(10) << "t" << std::setw(14) << "Euler"
              << std::setw(14) << "Exact" << std::setw(14) << "Error" << "\n";

    for (std::size_t i = 0; i < ys.size(); i += 100) {
        double t = t0 + i * h;
        double exact = y0 * std::exp(-t);
        double error = std::abs(ys[i] - exact);
        std::cout << std::left << std::setw(10) << t << std::setw(14) << ys[i]
                  << std::setw(14) << exact << std::setw(14) << error << "\n";
    }

    return 0;
}