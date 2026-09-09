#include <cassert>
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
    const double xi = 0.0;
    const double rho = -0.7;
    const double T = 1.0;
    const double h = 0.001;
    const int N = 50000;

    Eigen::Vector2d y0 = Eigen::Vector2d(S0, v0);

    auto driftS = [&](double t, Eigen::Vector2d y) { return r * y(0); };
    auto diffS = [&](double t, Eigen::Vector2d y) { return std::sqrt(y(1)) * y(0); };
    auto driftV = [&](double t, Eigen::Vector2d y) { return kappa * (theta - y(1)); };
    auto diffV = [&](double t, Eigen::Vector2d y) { return xi * std::sqrt(y(1)); };

    // Check the Feller condition
    if (2*kappa*theta < xi*xi) {
        std::cerr << "Warning: Feller condition violated\n";
    }

    double sum = 0;
    double sumSqr = 0;
    for (int i = 0; i < N; i++) {
        odelib::HestonEulerMaruyamaStepper<double> stepper(rho, i);
        auto last = odelib::integrate(stepper, driftS, diffS, driftV, diffV, y0, 0.0, T, h).back();
        sum += last(0);
        sumSqr += last(0) * last(0);
    }
    double mean = sum / N;

    // Var(X) = E[X^2] - (E[X])^2
    double singleVariance = sumSqr / N - mean * mean;

    // Var(X_bar) where X_bar = (1/N)(X1+...+XN):
    //   Var(X1+...+XN) = N * Var(X)      [variances add, since paths are independent]
    //   Var((1/N)*sum) = (1/N)^2 * N*Var(X) = Var(X)/N
    double meanVariance = singleVariance / N;

    double exact = S0 * exp(r * T);
    double error = std::abs(mean - exact);

    std::cout << "Average S_T: " << mean << "\n";
    std::cout << "Expected:    " << exact << "\n";
    std::cout << "Error:       " << error << "\n";

    assert(error < 2 * std::sqrt(meanVariance) && "Heston (xi=0) average S_T deviates too far from GBM expectation");
}