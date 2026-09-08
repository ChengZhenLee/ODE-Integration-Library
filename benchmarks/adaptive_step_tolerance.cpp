#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>

#include "odelib/integrate.hpp"
#include "odelib/steppers/adaptive_rk4.hpp"
#include "odelib/steppers/adaptive_rkdp.hpp"


struct BenchmarkPoint {
    public:
        double totalSteps;
        double error;
        double timeSeconds;
        double tolerance;

        BenchmarkPoint(double totalSteps, double error, double timeSeconds, double tolerance) :
            totalSteps(totalSteps), error(error), timeSeconds(timeSeconds), tolerance(tolerance) {}
};

auto steep_decay = [](double t, double y) { return -y * 50; };

int main(void) {
    const int repeats = 1000;
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;
    std::vector<double> tolerances = { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };

    double exact = y0 * std::exp(-t1 * 50);

    auto adaptive_RK4 = [&]() -> std::vector<BenchmarkPoint> {
        std::vector<BenchmarkPoint> points;

        for (const auto tolerance : tolerances) {
            double totalSteps;
            double result;
            double error;
            double duration;
            odelib::AdaptiveRK4Stepper<double, double> stepper(tolerance);

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < repeats; i++) {
                auto adaptiveIntegrationResult = odelib::integrate(stepper, steep_decay, y0, t0, t1, h);

                result = adaptiveIntegrationResult.states.back();
                totalSteps = adaptiveIntegrationResult.hs.size();
            }
            auto end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration<double>(end - start).count();

            error = std::abs(exact - result);
            points.push_back(BenchmarkPoint(totalSteps, error, duration, tolerance));
        }

        return points;
    };

    auto adaptive_RKDP = [&]() -> std::vector<BenchmarkPoint> {
        std::vector<BenchmarkPoint> points;

        for (const auto tolerance : tolerances) {
            double totalSteps;
            double result;
            double error;
            double duration;
            odelib::AdaptiveRKDPStepper<double, double> stepper(tolerance);

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < repeats; i++) {
                auto adaptiveIntegrationResult = odelib::integrate(stepper, steep_decay, y0, t0, t1, h);

                result = adaptiveIntegrationResult.states.back();
                totalSteps = adaptiveIntegrationResult.hs.size();
            }
            auto end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration<double>(end - start).count(); 

            error = std::abs(exact - result);
            points.push_back(BenchmarkPoint(totalSteps, error, duration, tolerance));
        }

        return points;
    };

    auto adaptiveRK4Future = std::async(std::launch::async, adaptive_RK4);
    auto adaptiveRKDPFuture = std::async(std::launch::async, adaptive_RKDP);

    auto adaptiveRK4Results = adaptiveRK4Future.get();
    auto adaptiveRKDPResults = adaptiveRKDPFuture.get();

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/data/AdaptiveStepTolerance.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }
    plotFile << "method,total_steps,error,time_seconds,tolerance" << "\n";
    
    for (const auto& point : adaptiveRK4Results) {
        plotFile << "AdaptiveRK4," << point.totalSteps << "," << point.error << "," << point.timeSeconds << "," << point.tolerance << "\n";
    }
    for (const auto& point : adaptiveRKDPResults) {
        plotFile << "AdaptiveRKDP," << point.totalSteps << "," << point.error << "," << point.timeSeconds << "," << point.tolerance << "\n";
    }

    plotFile.close();
}