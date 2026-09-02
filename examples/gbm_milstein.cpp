#include <cmath>
#include <iostream>
#include <iomanip>


#include "odelib/steppers/milstein.hpp"
#include "odelib/integrate.hpp"


double mu = 0.1;
double sigma = 0.2;

auto gbm_drift = [](double t, double y) { return mu * y; };
auto gbm_diffusion = [](double t, double y) { return sigma * y; };
auto gbm_diffusion_derivative = [](double t, double y) { return sigma; };

int main(void) {
    const double y0 = 1.0;   // initial condition
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;   // step size

    // Seed the stepper randomly
    odelib::MilsteinStepper<double> stepper(time(0));
    auto ys = odelib::integrate(stepper, gbm_drift, gbm_diffusion, gbm_diffusion_derivative, y0, t0, t1, h);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(10) << "t" << std::setw(16) << "Milstein" << "\n";

    for (std::size_t i = 0; i < ys.size(); i += 100) {
        double t = t0 + i * h;
        std::cout << std::left << std::setw(10) << t << std::setw(16) << ys[i] << "\n";
    }

    return 0;
}