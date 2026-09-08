#include <chrono>
#include <cmath>
#include <future>
#include <fstream>
#include <iostream>


#include "odelib/integrate.hpp"
#include "odelib/steppers/milstein.hpp"


struct BenchmarkPoint {
    int numThreads;
    double timeSeconds;
    double error;

    BenchmarkPoint(int numThreads, double timeSeconds, double error) :
        numThreads(numThreads), timeSeconds(timeSeconds), error(error) {}
};

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

    auto process_chunk = [&](int start, int end) -> double {
        double localSum = 0;
        double localSumSqr = 0;
        for (int i = start; i < end; i++) {
            odelib::MilsteinStepper<double, double> stepper(i);
            double S_T = odelib::integrate(stepper, drift, diffusion, diffusion_derivative, S0, 0.0, T, h).back();
            double payoff = std::max(S_T - K, 0.0);
            localSum += payoff;
        }
        return localSum;
    };

    std::vector<int> threadCount = { 1, 2, 4, 8 };
    std::vector<BenchmarkPoint> points;
    double exact = black_scholes(S0, K, r, sigma, T);
    
    for (auto numThreads : threadCount) {
        std::vector<std::future<double>> futures;

        double sum = 0;
        int chunkSize = N / numThreads;

        auto timeStart = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < numThreads; i++) {
            int start = i * chunkSize;
            int end = (i == numThreads - 1) ? N : i * chunkSize + chunkSize;
            futures.push_back(std::async(std::launch::async, process_chunk, start, end));
        }
        
        for (auto& future : futures) {
            sum += future.get();
        }

        double result = std::exp(-r * T) * sum / N;
        double error = std::abs(exact - result);

        auto timeEnd = std::chrono::high_resolution_clock::now();
        double timeSeconds = std::chrono::duration<double>(timeEnd - timeStart).count();

        points.push_back(BenchmarkPoint(numThreads, timeSeconds, error));
    }

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/data/ThreadSpeedup.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }

    plotFile << "method,thread_count,time_seconds,error\n";
    for (auto& point : points) {
        plotFile << "Milstein," << point.numThreads << "," << point.timeSeconds << "," << point.error << "\n";
    }
    plotFile.close();
}