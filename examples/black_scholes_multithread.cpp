#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>
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
    const int N = 100000;

    auto drift = [r](double t, double S) { return r * S; };
    auto diffusion = [sigma](double t, double S) { return sigma * S; };
    auto diffusion_derivative = [sigma](double t, double S) { return sigma; };

    auto processChunk = [&](int start, int end) -> std::pair<double, double> {
        double localSum = 0;
        double localSumSqr = 0;
        for (int i = start; i < end; i++) {
            odelib::MilsteinStepper<double, double> stepper(i);
            double S_T = odelib::integrate(stepper, drift, diffusion, diffusion_derivative, S0, 0.0, T, h).back();
            double payoff = std::max(S_T - K, 0.0);
            localSum += payoff;
            localSumSqr += payoff * payoff;
        }
        return {localSum, localSumSqr};
    };

    std::vector<std::future<std::pair<double, double>>> futures;
    int numThreads = std::thread::hardware_concurrency();
    int chunkSize = N / numThreads;
    
    for (int i = 0; i < numThreads; i++) {
        int start = i * chunkSize;
        int end = (i == numThreads - 1) ? N : start + chunkSize; // Last thread takes all remainder

        // Use std::launch::async to force each task onto its own thread immediately (not deffered)
        // also relies on having enough CPU cores available
        futures.push_back(std::async(std::launch::async, processChunk, start, end));
    }

    double sum = 0;
    double sumSqr = 0;

    for (auto& result : futures) {
        auto [localSum, localSumSqr] = result.get();
        sum += localSum;
        sumSqr += localSumSqr;
    }

    double price = std::exp(-r*T) * (sum / N);
    double variance = (std::exp(-r*T) * std::exp(-r*T)) * (sumSqr/N - sum/N * sum/N);
    double standardDeviation = std::sqrt(variance / N);

    double exact = black_scholes(S0, K, r, sigma, T);
    double error = std::abs(price - exact);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(16) << "Monte Carlo" << std::setw(16) << "Black Scholes"
              << std::setw(12) << "Error" << std::setw(12) << "Tolerance" << "\n";
    std::cout << std::left << std::setw(16) << price << std::setw(16) << exact
              << std::setw(12) << error << std::setw(12) << standardDeviation << "\n";
    std::cout << (error < standardDeviation ? "Match\n" : "Fail\n");

    return 0;
}