#include <fstream>
#include <future>
#include <iostream>

#include "odelib/integrate.hpp"
#include "odelib/steppers/adaptive_rk4.hpp"
#include "odelib/steppers/adaptive_rkdp.hpp"


auto steep_decay = [](double t, double y) { return -y * 50; };

int main(void) {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 5.0;
    const double h = 0.01;

    auto adaptive_RK4 = [&]() -> odelib::AdaptiveIntegrationResult<double, double> {
        odelib::AdaptiveRK4Stepper<double, double> stepper;
        return odelib::integrate(stepper, steep_decay, y0, t0, t1, h);
    };

    auto adaptive_RKDP = [&]() -> odelib::AdaptiveIntegrationResult<double, double> {
        odelib::AdaptiveRKDPStepper<double, double> stepper;
        return odelib::integrate(stepper, steep_decay, y0, t0, t1, h);
    };

    auto adaptiveRK4Future = std::async(std::launch::async, adaptive_RK4);
    auto adaptiveRKDPFuture = std::async(std::launch::async, adaptive_RKDP);

    auto adaptiveRK4Results = adaptiveRK4Future.get();
    auto adaptiveRKDPResults = adaptiveRKDPFuture.get();

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/data/AdaptiveStepTrace.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }
    plotFile << "method,t,h" << "\n";
    
    for (int i = 0; i < adaptiveRK4Results.states.size(); i++) {
        plotFile << "AdaptiveRK4," << adaptiveRK4Results.times[i] << "," << adaptiveRK4Results.hs[i] << "\n";
    }
    for (int i = 0; i < adaptiveRKDPResults.states.size(); i++) {
        plotFile << "AdaptiveRKDP," << adaptiveRKDPResults.times[i] << "," << adaptiveRKDPResults.hs[i] << "\n";
    }

    plotFile.close();
}