#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/adaptive_rk4.hpp"
#include "odelib/integrate.hpp"
#include "Eigen/Dense"


auto harmonic_oscillator = [](double t, Eigen::Vector2d y) -> Eigen::Vector2d { return Eigen::Vector2d(y(1), -y(0)); };

int main() {
    // y = (x, v)
    const Eigen::Vector2d y0(1.0, 1.0);
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;

    odelib::AdaptiveRK4Stepper<Eigen::Vector2d, double> stepper;
    auto result = odelib::integrate(stepper, harmonic_oscillator, y0, t0, t1, h);
    const auto& ys = result.states;
    const auto& ts = result.times;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "t\tAdaptive RK4\t\tExact\t\tError\n";

    for (std::size_t i = 0; i < ys.size(); i += 50) {
        double t = ts[i];
        Eigen::Vector2d exact = Eigen::Vector2d(y0(0) * cos(t) + y0(1) * sin(t), -y0(0) * sin(t) + y0(1) * cos(t));
        Eigen::Vector2d error = (ys[i] - exact).cwiseAbs();
        std::cout << t << "\t" << ys[i].transpose() << "\t" << exact.transpose() << "\t" << error.transpose()
                  << "\n";
    }

    return 0;
}