#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <random>

#include "odelib/integrate.hpp"
#include "odelib/steppers/milstein.hpp"
#include "odelib/steppers/euler_maruyama.hpp"
#include "test_utils.hpp"


double mu = 0.1;
double sigma = 0.2;

auto gbm_drift = [](double t, double y) { return mu * y; };
auto gbm_diffusion = [](double t, double y) { return sigma * y; };
auto gbm_diffusion_derivative = [](double t, double y) { return sigma; };

int main() {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 1.0;
    const double h = 0.001;
    const int N = 10000;

    TestUtils::printHeader("Monte Carlo Average Test");

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        odelib::MilsteinStepper<double> stepper(i);
        sum += odelib::integrate(stepper, gbm_drift, gbm_diffusion, gbm_diffusion_derivative, y0, t0, t1, h).back();
    }

    double expected = y0 * std::exp(mu * t1);
    double average = sum / static_cast<double>(N);
    double mcError = std::abs(expected - average);

    TestUtils::printRow("Calculated average:", average);
    TestUtils::printRow("Expected value:", expected);
    TestUtils::printRow("Error:", mcError);

    assert(mcError < 1.0 / std::sqrt(N) && "The average value deviated too far from the expected value");
    TestUtils::printPassed("Monte Carlo Average Test");


    TestUtils::printHeader("Strong Convergence Test");

    // Testing strong convergence
    double totalEMError = 0.0;
    double totalMError = 0.0;

    // We test for N number of independent specific paths
    for (int i = 0; i < N; i++) {
        std::mt19937_64 engine(i);
        std::normal_distribution<double> distribution(0, std::sqrt(h));

        // Passing 'i' into the constructor here is to prevent the 'most vexing parse'
        // dW is not generated in the Stepper, but passed in explicitly to ensure 'equality'
        odelib::EulerMaruyamaStepper<double> stepperEM(i);
        odelib::MilsteinStepper<double> stepperM(i);
        double yEM = y0;
        double yM = y0;
        double t = t0;
        double W = 0;

        while (t < t1) {
            double dW = distribution(engine);

            // Feed each step function with the same generated dW
            yEM = stepperEM.step(gbm_drift, gbm_diffusion, yEM, t, h, dW);
            yM = stepperM.step(gbm_drift, gbm_diffusion, gbm_diffusion_derivative, yM, t, h, dW);
            W += dW;
            t += h;
        }

        double exact = y0 * std::exp((mu - sigma*sigma/2)*t1 + sigma*W);
        double errorEM = std::abs(yEM - exact);
        double errorM = std::abs(yM - exact);

        totalEMError += errorEM;
        totalMError += errorM;
    }

    double avgEMError = totalEMError / static_cast<double>(N);
    double avgMError = totalMError / static_cast<double>(N);

    TestUtils::printRow("Avg EM error:", avgEMError);
    TestUtils::printRow("Avg Milstein error:", avgMError);

    assert(totalEMError > totalMError && "The error of strong convergence for the Milstein scheme should be smaller than Euler-Maruyama");
    TestUtils::printPassed("Strong Convergence Test");

    std::cout << "\nAll test_gbm_milstein checks passed.\n";
    return 0;
}
