#include <fstream>
#include <future>
#include <iostream>
#include <random>

#include "odelib/integrate.hpp"
#include "odelib/steppers/euler_maruyama.hpp"
#include "odelib/steppers/milstein.hpp"


double mu = 0.1;
double sigma = 0.2;

struct BenchmarkPoint {
    double h;
    double localSumEM;
    double localSumMilstein;

    BenchmarkPoint(double h, double localSumEM, double localSumMilstein) :
        h(h), localSumEM(localSumEM), localSumMilstein(localSumMilstein) {}
};

auto gbm_drift = [](double t, double y) { return mu * y; };
auto gbm_diffusion = [](double t, double y) { return sigma * y; };
auto gbm_diffusion_derivative = [](double t, double y) { return sigma; };


int main(void) {
    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t1 = 5.0;
    const std::vector<double> hs = {0.5, 0.1, 0.05, 0.01, 0.005, 0.001};
    const int N = 10000;

    auto process_chunk = [&](int start, int end) -> std::vector<BenchmarkPoint> {
        std::vector<BenchmarkPoint> res;
        for (double h : hs) {
            double localSumEM = 0;
            double localSumMilstein = 0;
            for (int i = start; i < end; i++) {
                // Make sure that the exact solution and steppers use the same underlying engine
                std::mt19937_64 engine(i);
                std::normal_distribution<double> distribution(0, std::sqrt(h));

                double y_em = y0, y_milstein = y0;
                double t = t0;
                double W = 0;
                odelib::EulerMaruyamaStepper<double, double> em_stepper;
                odelib::MilsteinStepper<double, double> milstein_stepper;

                while (t < t1) {
                    double dW = distribution(engine);
                    y_em = em_stepper.step(gbm_drift, gbm_diffusion, y_em, t, h, dW);
                    y_milstein = milstein_stepper.step(gbm_drift, gbm_diffusion, gbm_diffusion_derivative, y_milstein, t, h, dW);
                    W += dW;
                    t += h;
                }

                double exact = y0 * std::exp((mu - sigma*sigma/2)*t1 + sigma*W);
                localSumEM += std::abs(exact - y_em);
                localSumMilstein += std::abs(exact - y_milstein);
            }
            res.push_back(BenchmarkPoint(h, localSumEM, localSumMilstein));
        }
        return res;
    };

    int numThreads = std::thread::hardware_concurrency();
    int chunkSize = N / numThreads;
    std::vector<std::future<std::vector<BenchmarkPoint>>> futures;
    std::vector<double> eulerMaruyamaSums(hs.size());
    std::vector<double> milsteinSums(hs.size());

    for (int i = 0; i < numThreads; i++) {
        int start = i * chunkSize;
        int end = (i == numThreads - 1) ? N : i * chunkSize + chunkSize;
        futures.push_back(std::async(std::launch::async, process_chunk, start, end));
    }

    for (auto& future : futures) {
        auto points = future.get();
        for (size_t i = 0; i < hs.size(); i++) {
            eulerMaruyamaSums[i] += points[i].localSumEM;
            milsteinSums[i] += points[i].localSumMilstein;
        }
    }

    std::string outputPath = std::string(PROJECT_ROOT_DIR) + "/benchmarks/plots/StrongConvergencePlot.csv";
    std::ofstream plotFile(outputPath);
    if (!plotFile.is_open()) {
        std::cerr << "Failed to open output file!\n";
        return 1;
    }
    plotFile << "method,h,error\n";
    for (size_t i = 0; i < hs.size(); i++) {
        plotFile << "Euler-Maruyama," << hs[i] << "," << eulerMaruyamaSums[i]/N << "\n";
        plotFile << "Milstein," << hs[i] << "," << milsteinSums[i]/N << "\n";
    }
    plotFile.close();
}
