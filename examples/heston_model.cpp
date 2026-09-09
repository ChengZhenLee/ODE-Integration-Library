#include <iomanip>
#include <iostream>


#include "odelib/integrate.hpp"
#include "odelib/steppers/heston_euler_maruyama.hpp"
#include "Eigen/Dense"


int main(void) {
    const double S0 = 100.0;
    const double v0 = 0.04;
    const double r = 0.05;
    const double kappa = 2.0;
    const double theta = 0.04;
    const double xi = 0.3;
    const double rho = -0.7;
    const double T = 1.0;
    const double h = 0.001;

    Eigen::Vector2d y0 = Eigen::Vector2d(S0, v0);

    auto driftS = [&](double t, Eigen::Vector2d y) { return r * y(0); };
    auto diffS = [&](double t, Eigen::Vector2d y) { return std::sqrt(y(1)) * y(0); };
    auto driftV = [&](double t, Eigen::Vector2d y) { return kappa * (theta - y(1)); };
    auto diffV = [&](double t, Eigen::Vector2d y) { return xi * std::sqrt(y(1)); };

    // Check the Feller condition
    if (2*kappa*theta < xi*xi) {
        std::cerr << "Warning: Feller condition violated\n";
    }

    odelib::HestonEulerMaruyamaStepper<double> stepper(rho, 0);
    auto result = odelib::integrate(stepper, driftS, diffS, driftV, diffV, y0, 0.0, T, h);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "t\t\tS\t\tv\n";

    // print every 100th point so the table stays a manageable size
    for (std::size_t i = 0; i < result.size(); i += 100) {
        double t = static_cast<double>(i) * h;
        std::cout << t << "\t\t" << result[i](0) << "\t\t" << result[i](1) << "\n";
    }

    // always print the last point
    std::size_t last = result.size() - 1;
    double t_last = static_cast<double>(last) * h;
    std::cout << t_last << "\t\t" << result[last](0) << "\t\t" << result[last](1) << "\n";
}