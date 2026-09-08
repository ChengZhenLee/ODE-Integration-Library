#include <cmath>
#include <fstream>
#include <iostream>
#include <future>


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
    std::vector<int> Ns = {1000, 5000, 10000, 50000, 100000, 500000};

    auto drift = [r](double t, double S) { return r * S; };
    auto diffusion = [sigma](double t, double S) { return sigma * S; };
    auto diffusion_derivative = [sigma](double t, double S) { return sigma; };

    auto processChunk = [&](int n) -> double {
        double sum = 0;
        for (int i = 0; i < n; i++) {
            odelib::MilsteinStepper<double, double> stepper(i);
            double result = odelib::integrate(stepper, drift, diffusion, diffusion_derivative, S0, 0.0, T, h).back();
            sum += std::max(result - K, 0.0);
        }
        return exp(-r*T) * sum/n;
    };
    
    std::vector<std::future<double>> futures;
    std::vector<double> errors;

    for (int N : Ns) {
        futures.push_back(std::async(std::launch::async, processChunk, N));
    }

    double exact = black_scholes(S0, K, r, sigma, T);
    for (size_t i = 0; i < Ns.size(); i++) {
        errors.push_back(std::abs(exact - futures[i].get()));
    }

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/plots/MonteCarloN.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }

    plotFile << "method,N,error\n";
    for (size_t i = 0; i < Ns.size(); i++) {
        plotFile << "Milstein," << Ns[i] << "," << errors[i] << "\n";
    }
    plotFile.close();
}