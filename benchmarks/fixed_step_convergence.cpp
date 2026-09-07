#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <future>

#include "odelib/integrate.hpp"
#include "odelib/steppers/euler.hpp"
#include "odelib/steppers/rk4.hpp"
#include "odelib/steppers/rkdp.hpp"

struct BenchmarkPoint {
    public:
        double h;
        double error;
        double timeSeconds;

        BenchmarkPoint(double h, double error, double timeSeconds) : 
            h(h), error(error), timeSeconds(timeSeconds) {}
};

auto decay = [](double t, double y) { return -y; };

int main(void) {
    const int repeats = 100;
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 5.0;
    const std::vector<double> hs = {0.5, 0.1, 0.05, 0.01, 0.005, 0.001};

    double exact = y0 * std::exp(-t1);

    auto euler = [&]() -> std::vector<BenchmarkPoint> {
        std::vector<BenchmarkPoint> res;
        for (const auto h : hs) {
            double result;
            double bestTime = std::numeric_limits<double>::max();
            double error;
            odelib::EulerStepper<double, double> stepper;

            for (int i = 0; i < repeats; i++) {
                auto start = std::chrono::high_resolution_clock::now();
                result = odelib::integrate(stepper, decay, y0, t0, t1, h).back();
                auto end = std::chrono::high_resolution_clock::now();
                bestTime = std::min(bestTime, std::chrono::duration<double>(end - start).count());
            }

            error = std::abs(result - exact);
            res.push_back(BenchmarkPoint(h, error, bestTime));
        }

        return res;
    };

    auto rk4 = [&]() -> std::vector<BenchmarkPoint> {
        std::vector<BenchmarkPoint> res;
        for (const auto h : hs) {
            double result;
            double bestTime = std::numeric_limits<double>::max();
            double error;
            odelib::RK4Stepper<double, double> stepper;

            for (int i = 0; i < repeats; i++) {
                auto start = std::chrono::high_resolution_clock::now();
                result = odelib::integrate(stepper, decay, y0, t0, t1, h).back();
                auto end = std::chrono::high_resolution_clock::now();
                bestTime = std::min(bestTime, std::chrono::duration<double>(end - start).count());
            }

            error = std::abs(result - exact);
            res.push_back(BenchmarkPoint(h, error, bestTime));
        }

        return res;
    };

    auto rkdp = [&]() -> std::vector<BenchmarkPoint> {
        std::vector<BenchmarkPoint> res;
        for (const auto h : hs) {
            double result;
            double bestTime = std::numeric_limits<double>::max();
            double error;
            odelib::RKDPStepper<double, double> stepper;

            for (int i = 0; i < repeats; i++) {
                auto start = std::chrono::high_resolution_clock::now();
                result = odelib::integrate(stepper, decay, y0, t0, t1, h).back();
                auto end = std::chrono::high_resolution_clock::now();
                bestTime = std::min(bestTime, std::chrono::duration<double>(end - start).count());
            }

            error = std::abs(result - exact);
            res.push_back(BenchmarkPoint(h, error, bestTime));
        }

        return res;
    };

    auto futureEuler = std::async(std::launch::async, euler);
    auto futureRK4 = std::async(std::launch::async, rk4);
    auto futureRKDP = std::async(std::launch::async, rkdp);

    auto eulerResults = futureEuler.get();
    auto rk4Results = futureRK4.get();
    auto rkdpResults = futureRKDP.get();

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/plots/FixedStepConvergence.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }
    plotFile << "method,h,error,time_seconds\n";

    for (auto& result : eulerResults) {
        plotFile << "Euler," << result.h << "," << result.error << "," << result.timeSeconds << "\n";
    }
    for (auto& result : rk4Results) {
        plotFile << "RK4," << result.h << "," << result.error << "," << result.timeSeconds << "\n";
    }
    for (auto& result : rkdpResults) {
        plotFile << "RKDP," << result.h << "," << result.error << "," << result.timeSeconds << "\n";
    }

    plotFile.close();
}