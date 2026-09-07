#include <cmath>
#include <iostream>


#include "odelib/integrate.hpp"
#include "odelib/steppers/milstein.hpp"


double phi(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double black_scholes(double S0, double K, double r, double sigma, double T) {
    double d1 = (std::log(S0/K) + (r + sigma*sigma/2)*T) / (sigma*std::sqrt(T));
    double d2 = d1 - sigma*std::sqrt(T);
    return S0*phi(d1) - K*std::exp(-r*T)*phi(d2);
}



int main(void) {
    const double S0 = 100;
    const double K = 100;
    const double r = 0.05;
    const double sigma = 0.2;
    const double T = 1.0;
    const double h = 0.01;
    const int N = 100000;

    auto drift = [r](double t, double S) { return r * S; };
    auto diffusion = [sigma](double t, double S) { return sigma * S; };
    auto diffusion_derivative = [sigma](double t, double S) { return sigma; };

    double sum = 0;
    double sumSqr = 0;
    for (int i = 0; i < N; i++) {
        odelib::MilsteinStepper<double, double> stepper(i);
        double S_T = odelib::integrate(stepper, drift, diffusion, diffusion_derivative, S0, 0.0, T, h).back();
        sum += std::max(S_T - K, 0.0);
        sumSqr += sum * sum;
    }
    double price = std::exp(-r*T) * (sum / N);
    double variance = (std::exp(-r*T) * std::exp(-r*T)) * (sumSqr/N - sum/N * sum/N);
    double standardDeviation = std::sqrt(variance / N);

    double exact = black_scholes(S0, K, r, sigma, T);
    double error = std::abs(price - exact);

    std::cout << price << "        " << exact << "        " << error << "\n";
    if (error < standardDeviation/std::sqrt(N)) {
        std::cout << "Match\n";
    } else {
        std::cout << "Fail\n";
    }

    return 0;
}