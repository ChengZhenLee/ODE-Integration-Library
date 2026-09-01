#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/adaptive_rk4.hpp"
#include "odelib/integrate.hpp"


double k = 0.5;
auto scaled_decay = [](double t, double y) { return -k * y; };

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;

    odelib::AdaptiveRK4Stepper<double, double> stepper;
    auto result = odelib::integrate(stepper, scaled_decay, y0, t0, t1, h);
    const auto& ys = result.states;
    const auto& ts = result.times;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "t\tAdaptive RK4\t\tExact\t\tError\n";

    for (std::size_t i = 0; i < ys.size(); i ++) {
        double t = ts[i];
        double exact = y0 * std::exp(-t * k);
        double error = std::abs(ys[i] - exact);
        std::cout << t << "\t" << ys[i] << "\t" << exact << "\t" << error
                  << "\n";
    }

    return 0;
}