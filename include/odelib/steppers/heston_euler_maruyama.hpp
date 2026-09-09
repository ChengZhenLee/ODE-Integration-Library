#pragma once

#include <random>
#include "Eigen/Dense"


namespace odelib {

template <typename Time = double>
class HestonEulerMaruyamaStepper {
    private:
        std::mt19937_64 engine;
        double rho = -0.4;

    public:
        HestonEulerMaruyamaStepper() {}

        HestonEulerMaruyamaStepper(double rho, unsigned int seed) : rho(rho), engine(seed) {}

        template <typename DriftS, typename DiffS, typename DriftV, typename DiffV>
        Eigen::Vector2d step(DriftS&& driftS, DiffS&& diffS, DriftV&& driftV, DiffV&& diffV, 
            const Eigen::Vector2d& y, Time t, Time h) {
            std::normal_distribution<double> distribution(0, 1);
            
            double dW1 = std::sqrt(h) * distribution(engine);
            double dW2 = std::sqrt(h) * (rho * distribution(engine) + std::sqrt(1 - rho*rho) * distribution(engine));

            double sNext = y(0) + driftS(t, y) * h + diffS(t, y) * dW1;

            // Ensure the variance is never negative
            double vNext = std::max(y(1) + driftV(t, y) * h + diffV(t, y) * dW2, 0.0);
            return Eigen::Vector2d(sNext, vNext);
        }
};
}