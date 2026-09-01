#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/rk4.hpp"
#include "odelib/integrate.hpp"
#include "Eigen/Dense"


auto decay = [](double t, Eigen::Vector2d y) { return -y; };

int main() {
    const Eigen::Vector2d y0(1.0, 2.0);
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;

    odelib::RK4Stepper<Eigen::Vector2d, double> stepper;
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "t\tRK4\t\tExact\t\tError\n";

    for (std::size_t i = 0; i < ys.size(); i += 100) {
        double t = t0 + i * h;
        Eigen::Vector2d exact = y0 * std::exp(-t);
        Eigen::Vector2d error = (ys[i] - exact).cwiseAbs();
        std::cout << t << "\t" << ys[i].transpose() << "\t" << exact.transpose() << "\t" << error.transpose()
                  << "\n";
    }

    return 0;
}