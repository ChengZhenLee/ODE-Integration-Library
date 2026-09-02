#pragma once


#include <random>


namespace odelib {

template <typename State, typename Time=double>
struct MilsteinStepper {
    private:
        std::mt19937_64 engine;

    public:
        MilsteinStepper(unsigned int seed) : engine(seed) {}

        template <typename Drift, typename Diffusion, typename DiffusionDerivative>
        State step(Drift&& f, Diffusion&& g, DiffusionDerivative&& gPrime, State y, Time t, Time h) {
            std::normal_distribution<Time> distribution(0, std::sqrt(h));
            Time dW = distribution(engine);

            return y + f(t, y) * h + g(t, y) * dW + 0.5 * g(t, y) * gPrime(t, y) * (std::pow(dW, 2) - h);
        }

        template <typename Drift, typename Diffusion, typename DiffusionDerivative>
        State step(Drift&& f, Diffusion&& g, DiffusionDerivative&& gPrime, State y, Time t, Time h, Time dW) {
            return y + f(t, y) * h + g(t, y) * dW + 0.5 * g(t, y) * gPrime(t, y) * (std::pow(dW, 2) - h);
        }
};

}