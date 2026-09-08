#include <cmath>
#include <fstream>
#include <iostream>
#include <future>
#include <random>


#include "odelib/integrate.hpp"
#include "odelib/steppers/euler_maruyama.hpp"
#include "odelib/steppers/milstein.hpp"


struct BenchmarkPoint {
    double t;
    double yEM;
    double yMilstein;
    double yExact;

    BenchmarkPoint(double t, double yEM, double yMilstein, double yExact) : 
        t(t), yEM(yEM), yMilstein(yMilstein), yExact(yExact) {}
};


int main(void) {
    const double S0 = 100;
    const double K = 100;
    const double r = 0.05;
    const double sigma = 0.2;
    const double T = 10.0;
    const double h = 0.5;

    auto drift = [r](double t, double S) { return r * S; };
    auto diffusion = [sigma](double t, double S) { return sigma * S; };
    auto diffusion_derivative = [sigma](double t, double S) { return sigma; };

    double t0 = 0;
    double W = 0;
    double yEM = S0;
    double yMilstein = S0;
    std::mt19937_64 engine(0);
    std::normal_distribution<double> distribution(0, std::sqrt(h));
    odelib::EulerMaruyamaStepper<double, double> stepperEM;
    odelib::MilsteinStepper<double, double> stepperMilstein;
    std::vector<BenchmarkPoint> points = {BenchmarkPoint(t0, yEM, yMilstein, S0)};

    long n = static_cast<long>(std::round((T - t0) / h));
    for (long i = 0; i < n; i++) {
        double t = t0 + static_cast<double>(i) * h;
        double dW = distribution(engine);
        yEM = stepperEM.step(drift, diffusion, yEM, t, h, dW);
        yMilstein = stepperMilstein.step(drift, diffusion, diffusion_derivative, yMilstein, t, h, dW);
        W += dW;
        double tNext = t + h;
        double exact = S0 * std::exp((r -  sigma * sigma / 2) * tNext + sigma * W);
        points.push_back(BenchmarkPoint(tNext, yEM, yMilstein, exact));
    };

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/data/WalkComparison.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }

    plotFile << "method,t,state\n";
    for (auto& point : points) {
        plotFile << "Euler-Maruyama," << point.t << "," << point.yEM << "\n";
        plotFile << "Milstein," << point.t << "," << point.yMilstein << "\n";
        plotFile << "Exact," << point.t << "," << point.yExact << "\n";
    }
    plotFile.close();
}