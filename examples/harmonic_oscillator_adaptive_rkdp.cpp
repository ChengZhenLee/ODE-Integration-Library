#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/adaptive_rkdp.hpp"
#include "odelib/integrate.hpp"
#include "Eigen/Dense"


auto harmonic_oscillator = [](double t, Eigen::Vector2d y) -> Eigen::Vector2d { return Eigen::Vector2d(y(1), -y(0)); };

int main() {
    // y = (x, v)
    const Eigen::Vector2d y0(1.0, 1.0);
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;

    odelib::AdaptiveRKDPStepper<Eigen::Vector2d, double> stepper;
    auto result = odelib::integrate(stepper, harmonic_oscillator, y0, t0, t1, h);
    const auto& ys = result.states;
    const auto& ts = result.times;
    const auto& hs = result.hs;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "t\tAdaptive RKDP\t\tStep Size\t\tExact\t\tError\n";

    for (std::size_t i = 0; i < ys.size(); i += 10) {
        double t = ts[i];
        Eigen::Vector2d exact = Eigen::Vector2d(y0(0) * cos(t) + y0(1) * sin(t), -y0(0) * sin(t) + y0(1) * cos(t));
        Eigen::Vector2d error = (ys[i] - exact).cwiseAbs();
        std::cout << t << "\t" << ys[i].transpose() << "\t" << hs[i] << "\t" << exact.transpose() << "\t" << error.transpose()
                  << "\n";
    }

    return 0;
}