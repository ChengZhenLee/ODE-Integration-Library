#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "odelib/integrate.hpp"
#include "odelib/steppers/euler_maruyama.hpp"
#include "test_utils.hpp"


double mu = 0.1;
double sigma = 0.2;

auto gbm_drift = [](double t, double y) { return mu * y; };
auto gbm_diffusion = [](double t, double y) { return sigma * y; };
auto gbm_diffusion_empty = [](double t, double y) { return 0 * y; };

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.001;
    const int N = 10000;

    TestUtils::printHeader("Monte Carlo Average Test");

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        odelib::EulerMaruyamaStepper<double> stepper(i);
        sum += odelib::integrate(stepper, gbm_drift, gbm_diffusion, y0, t0, t1, h).back();
    }

    double expected = y0 * std::exp(mu * t1);
    double average = sum / static_cast<double>(N);
    double mcError = std::abs(expected - average);

    TestUtils::printRow("Calculated average:", average);
    TestUtils::printRow("Expected value:", expected);
    TestUtils::printRow("Error:", mcError);

    assert(mcError < 1.0 / std::sqrt(N) && "The average value deviated too far from the expected value");
    TestUtils::printPassed("Monte Carlo Average Test");

    TestUtils::printHeader("Deterministic Path Test");

    odelib::EulerMaruyamaStepper<double> stepper(static_cast<unsigned int>(time(0)));
    const auto& ys = odelib::integrate(stepper, gbm_drift, gbm_diffusion_empty, y0, t0, t1, h);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  " << std::left << std::setw(10) << "t" << std::setw(16) << "Euler-Maruyama"
              << std::setw(16) << "Exact" << std::setw(16) << "Error" << "\n";

    double maxError = 0.0;
    for (std::size_t i = 0; i < ys.size(); i += 50) {
        double t = t0 + i * h;
        double exact = std::exp(mu * t);
        double error = std::abs(exact - ys[i]);
        maxError = std::max(maxError, error);

        std::cout << "  " << std::left << std::setw(10) << t << std::setw(16) << ys[i]
                  << std::setw(16) << exact << std::setw(16) << error << "\n";
    }

    assert(maxError < 0.01 && "The non-probabilistic Euler-Maruyama solution is too far away from the exact solution");
    TestUtils::printPassed("Deterministic Path Test");

    std::cout << "\nAll test_gbm_euler_maruyama checks passed.\n";
    return 0;
}
