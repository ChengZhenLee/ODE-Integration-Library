#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/rk4.hpp"
#include "odelib/integrate.hpp"


double k = 0.5;
auto scaled_decay = [](double t, double y) { return -k * y; };

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;

    odelib::RK4Stepper<double> stepper;
    auto ys = odelib::integrate(stepper, scaled_decay, y0, t0, t1, h);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(10) << "t" << std::setw(14) << "RK4"
              << std::setw(14) << "Exact" << std::setw(14) << "Error" << "\n";

    for (std::size_t i = 0; i < ys.size(); i += 100) {
        double t = t0 + i * h;
        double exact = y0 * std::exp(-t * k);
        double error = std::abs(ys[i] - exact);
        std::cout << std::left << std::setw(10) << t << std::setw(14) << ys[i]
                  << std::setw(14) << exact << std::setw(14) << error << "\n";
    }

    return 0;
}