#include <cmath>
#include <fstream>
#include <iostream>
#include <future>
#include <random>


#include "odelib/integrate.hpp"
#include "odelib/steppers/euler_maruyama.hpp"
#include "odelib/steppers/milstein.hpp"


struct BenchmarkPoint {
    double yEM;
    double yMilstein;
    double yExact;

    BenchmarkPoint(double yEM, double yMilstein, double yExact) : 
        yEM(yEM), yMilstein(yMilstein), yExact(yExact) {}
};


int main(void) {
    const double S0 = 100;
    const double K = 100;
    const double r = 0.05;
    const double sigma = 0.2;
    const double T = 1.0;
    const double h = 0.01;

    auto drift = [r](double t, double S) { return r * S; };
    auto diffusion = [sigma](double t, double S) { return sigma * S; };
    auto diffusion_derivative = [sigma](double t, double S) { return sigma; };

    double t = 0;
    double W = 0;
    double yEM = S0;
    double yMilstein = S0;
    std::mt19937_64 engine(0);
    std::normal_distribution<double> distribution(0, std::sqrt(h));
    odelib::EulerMaruyamaStepper<double, double> stepperEM;
    odelib::MilsteinStepper<double, double> stepperMilstein;
    std::vector<BenchmarkPoint> points = {BenchmarkPoint(yEM, yMilstein, S0)};

    while (t < T) {
        double dW = distribution(engine);
        yEM = stepperEM.step(drift, diffusion, yEM, t, h, dW);
        yMilstein = stepperMilstein.step(drift, diffusion, diffusion_derivative, yMilstein, t, h, dW);
        W += dW;
        t += h;
        double exact = S0 * std::exp((r -  sigma * sigma / 2) * t + sigma * W);
        points.push_back(BenchmarkPoint(yEM, yMilstein, exact));
    };

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/plots/WalkComparison.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }

    plotFile << "method,state\n";
    for (auto& point : points) {
        plotFile << "Euler-Maruyama," << point.yEM << "\n";
        plotFile << "Milstein," << point.yMilstein << "\n";
        plotFile << "Exact," << point.yExact << "\n";
    }
    plotFile.close();
}