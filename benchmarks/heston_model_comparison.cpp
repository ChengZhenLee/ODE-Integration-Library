#include <fstream>
#include <iostream>
#include <vector>


#include "odelib/integrate.hpp"
#include "odelib/steppers/heston_euler_maruyama.hpp"
#include "Eigen/Dense"


int main(void) {
    const double S0 = 100.0;
    const double v0 = 0.04;
    const double r = 0.05;
    const double kappa = 2.0;
    const double theta = 0.04;
    const double xi = 0.3;
    const double rho = -0.7;
    const double T = 1.0;
    const double h = 0.001;

    const int repeats = 5;
    std::vector<std::vector<Eigen::Vector2d>> res;
    Eigen::Vector2d y0 = Eigen::Vector2d(S0, v0);

    auto driftS = [&](double t, Eigen::Vector2d y) { return r * y(0); };
    auto diffS = [&](double t, Eigen::Vector2d y) { return std::sqrt(y(1)) * y(0); };
    auto driftV = [&](double t, Eigen::Vector2d y) { return kappa * (theta - y(1)); };
    auto diffV = [&](double t, Eigen::Vector2d y) { return xi * std::sqrt(y(1)); };

    // Check the Feller condition
    if (2*kappa*theta < xi*xi) {
        std::cerr << "Warning: Feller condition violated\n";
    }

    for (int i = 0; i < repeats; i++) {
        odelib::HestonEulerMaruyamaStepper<double> stepper(rho, i);
        res.push_back(odelib::integrate(stepper, driftS, diffS, driftV, diffV, y0, 0.0, T, h));
    }

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/data/HestonModelComparison.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }

    plotFile << "method,walk,t,s,v\n";
    for (int i = 0; i < repeats; i++) {
        for (size_t j = 0; j < res[i].size(); j++) {
            plotFile << "Euler-Maruyama," << i << "," << static_cast<double>(j) * h << "," << res[i][j](0) << "," << res[i][j](1) << "\n";
        }
    }   
    plotFile.close();
}